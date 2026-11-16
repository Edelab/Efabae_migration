#To resolve the putative clusters of migrants within a month, 
#we adopt the same method as in Reich et al. 2021
#Here, we used the k-means to cluster those individuals, 
#we tested three transformations to see the best number of clusters
#This code was derived from Reich et al. 2021, with debug and some
#extensions assisted by the Claude Sonnet 5.5, including the organization
#to facilitate the script run.

#Reich et al. Methods Ecol Evol 12, 2445–2457 (2021).

library(terra)
library(ggplot2)
library(cluster)
library(readxl)

# -----------------------------------------------------------------------------
# 1. Settings
# -----------------------------------------------------------------------------
probability_dir     <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/cluster"   # parent folder; searched recursively
classification_xlsx <- "C:/Users/abraa/OneDrive/Documents/Empoasca/Isotope/2assign/cluster/local_migratory_classification.xlsx"
class_sheet         <- 1

id_col          <- "ID"
call_col        <- "Final_putative"
migratory_label <- "Migratory"
combo_col       <- "Combo"
er_col          <- "Enrichment_ratio"
mindist_col     <- "Min_distance"

output_dir <- file.path(probability_dir, "cluster_migrants_by_month")

months_to_run <- c("JUN", "JUL", "AUG", "SEP")
transforms    <- c("none", "sqrt", "rank")
agg_factor    <- 2        # we reduce the agregation to 2 to avoid memory issues
k_range       <- 2:4      # candidate k
min_n <- 4 


# NULL = best average silhouette; integer = override (the same k is used for all transforms)
k_by_month <- list(JUN = NULL, JUL = NULL, AUG = NULL, SEP = 2)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (tr in transforms) dir.create(file.path(output_dir, tr, "cluster_mean_surfaces"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Load samples
# -----------------------------------------------------------------------------
cls_all <- as.data.frame(read_excel(classification_xlsx, sheet = class_sheet))
miss <- setdiff(c(id_col, call_col), names(cls_all))
if (length(miss) > 0) stop("Column(s) not found in the classification file: ", paste(miss, collapse = ", "))
for (cc in c(combo_col, er_col, mindist_col)) if (!cc %in% names(cls_all)) cls_all[[cc]] <- NA

cat("Final call, all individuals:\n"); print(table(cls_all[[call_col]], useNA = "ifany"))
cls <- cls_all[!is.na(cls_all[[call_col]]) & cls_all[[call_col]] == migratory_label, , drop = FALSE]

id_pat <- "^A[0-9]+_([0-9]{1,2})([A-Z]{3})([0-9]{2})$"
cls <- cls[grepl(id_pat, cls[[id_col]]), , drop = FALSE]          # drops _DUP / non-A samples
cls$Month <- sub(id_pat, "\\2", cls[[id_col]])
cls$Year  <- paste0("20", sub(id_pat, "\\3", cls[[id_col]]))
cls$Combo <- as.character(cls[[combo_col]])
cls$ER    <- suppressWarnings(as.numeric(cls[[er_col]]))
cls$MinDist <- suppressWarnings(as.numeric(cls[[mindist_col]]))
cls$Has_H <- grepl("(^|_)H(_|$)", cls$Combo)

months <- intersect(months_to_run, unique(cls$Month))
cat(sprintf("\nMigratory individuals used: %d\n", nrow(cls)))
print(table(Month = cls$Month, Year = cls$Year))

# lookup of every probability tif under probability_dir: ID -> path
pf    <- list.files(probability_dir, pattern = "_probability\\.tif$", recursive = TRUE, full.names = TRUE)
pf_id <- sub("_probability\\.tif$", "", basename(pf))
lookup <- setNames(pf[!duplicated(pf_id)], pf_id[!duplicated(pf_id)])
cat(sprintf("Found %d probability tifs under %s\n", length(lookup), probability_dir))

# -----------------------------------------------------------------------------
# 3. Functions for adjustment in the bugs previously identified
# -----------------------------------------------------------------------------
make_mat <- function(V, tr) {                     # V = cells x individuals
  if (tr == "none") return(t(V))
  if (tr == "sqrt") return(t(sqrt(V)))
  t(apply(V, 2, function(p) rank(p, ties.method = "average") / length(p)))   # rank
}

ari <- function(a, b) {                           # adjusted Rand index
  tab <- table(a, b); n <- sum(tab)
  sc <- function(x) sum(choose(x, 2))
  sij <- sc(tab); sa <- sc(rowSums(tab)); sb <- sc(colSums(tab))
  ex <- sa * sb / choose(n, 2)
  if ((0.5 * (sa + sb) - ex) == 0) return(NA_real_)
  (sij - ex) / (0.5 * (sa + sb) - ex)
}

fmt_tab <- function(x) paste(names(table(x)), as.integer(table(x)), sep = ":", collapse = ", ")

# -----------------------------------------------------------------------------
# 4. Per-month clustering testing the three transformations
# -----------------------------------------------------------------------------
cluster_month <- function(m) {
  set.seed(1)
  sub_df <- cls[cls$Month == m, , drop = FALSE]
  files  <- unname(lookup[sub_df[[id_col]]])
  ok     <- !is.na(files)
  if (any(!ok)) cat(sprintf("  %s: no probability tif for: %s\n", m, paste(sub_df[[id_col]][!ok], collapse = ", ")))
  sub_df <- sub_df[ok, , drop = FALSE]; files <- files[ok]
  n <- nrow(sub_df)
  cat(sprintf("\n==== %s (n = %d) ====\n", m, n))
  if (n < min_n) { cat("  too few individuals, skipping.\n"); return(NULL) }
  
  # bring all surfaces onto one grid (resample any that differ)
  rl   <- lapply(files, rast)
  tmpl <- rl[[which.max(sapply(rl, function(x) prod(dim(x)[1:2])))]]
  bad  <- !sapply(rl, function(x) compareGeom(x, tmpl, stopOnError = FALSE))
  if (any(bad)) {
    cat(sprintf("  Note: different grid, resampled to template: %s\n", paste(sub_df[[id_col]][bad], collapse = ", ")))
    rl[bad] <- lapply(rl[bad], function(x) resample(x, tmpl, method = "bilinear"))
  }
  r <- rast(rl)
  if (agg_factor > 1) r <- aggregate(r, fact = agg_factor, fun = "mean", na.rm = TRUE)
  vals  <- values(r)
  valid <- complete.cases(vals)
  V     <- vals[valid, , drop = FALSE]
  V     <- sweep(V, 2, colSums(V), "/")             # each surface sums to 1
  rm(vals, rl)
  
  sub_df$Conc50 <- round(apply(V, 2, function(p) { p <- sort(p, decreasing = TRUE); which(cumsum(p) >= 0.5)[1] / length(p) }), 4)
  
  ks <- k_range[k_range <= max(2, floor(n / 2))]
  diag_list <- list(); summ_list <- list(); k_used <- c()
  
  for (tr in transforms) {
    cat(sprintf("\n-- %s | transform: %s --\n", m, tr))
    mat <- make_mat(V, tr); rownames(mat) <- sub_df[[id_col]]
    d <- dist(mat)
    
    wss <- numeric(length(ks)); sil <- numeric(length(ks))
    for (i in seq_along(ks)) {
      km_i   <- kmeans(mat, centers = ks[i], nstart = 25, iter.max = 50)
      wss[i] <- km_i$tot.withinss
      sil[i] <- mean(cluster::silhouette(km_i$cluster, d)[, 3])
    }
    print(data.frame(k = ks, avg_silhouette = round(sil, 3)), row.names = FALSE)
    diag_list[[tr]] <- rbind(
      data.frame(transform = tr, k = ks, metric = "Within-cluster SS (normalised)", value = wss / max(wss)),
      data.frame(transform = tr, k = ks, metric = "Average silhouette", value = sil))
    
    k_choice <- k_by_month[[m]]
    if (is.null(k_choice)) { k_choice <- ks[which.max(sil)]; cat(sprintf("  auto k = %d\n", k_choice)) }
    else cat(sprintf("  override k = %d\n", k_choice))
    km <- kmeans(mat, centers = k_choice, nstart = 25, iter.max = 50)
    sub_df[[paste0("Cluster_", tr)]] <- km$cluster
    k_used[tr] <- k_choice
    
    # mean surface per cluster (from the renormalised surfaces)
    for (cl in sort(unique(km$cluster))) {
      mv  <- rowMeans(V[, km$cluster == cl, drop = FALSE]); mv <- mv / sum(mv)
      out <- rast(r[[1]]); values(out) <- NA; out[which(valid)] <- mv
      writeRaster(out, file.path(output_dir, tr, "cluster_mean_surfaces", sprintf("%s_cluster%d_mean.tif", m, cl)), overwrite = TRUE)
    }
    
    # per-cluster summary
    cl_vec <- km$cluster
    summ_list[[tr]] <- do.call(rbind, lapply(sort(unique(cl_vec)), function(cl) {
      s <- sub_df[cl_vec == cl, ]
      data.frame(Month = m, Transform = tr, Cluster = cl, n = nrow(s),
                 Years = fmt_tab(s$Year),
                 Conc50_median = round(median(s$Conc50), 3),
                 ER_median = round(median(s$ER, na.rm = TRUE), 2),
                 ER_min = round(suppressWarnings(min(s$ER, na.rm = TRUE)), 2),
                 ER_max = round(suppressWarnings(max(s$ER, na.rm = TRUE)), 2),
                 MinDist_median_km = round(median(s$MinDist, na.rm = TRUE), 0),
                 Isotope_sets = fmt_tab(s$Combo), stringsAsFactors = FALSE)
    }))
    print(summ_list[[tr]][, c("Cluster", "n", "Years", "ER_median", "Isotope_sets")], row.names = FALSE)
    
    tb <- table(Cluster = cl_vec, HasH = sub_df$Has_H)
    if (ncol(tb) > 1) {
      ft <- tryCatch(fisher.test(tb), error = function(e) NULL)
      if (!is.null(ft)) cat(sprintf("  cluster vs has-hydrogen: Fisher p = %.3f\n", ft$p.value))
    }
  }
  
  # figure: k diagnostics, one panel per transform, no title
  dg <- do.call(rbind, diag_list); dg$transform <- factor(dg$transform, levels = transforms)
  p <- ggplot(dg, aes(x = k, y = value, color = metric)) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ transform, nrow = 1) +
    scale_x_continuous(breaks = ks) +
    scale_color_manual(values = c("Within-cluster SS (normalised)" = "steelblue", "Average silhouette" = "darkorange")) +
    labs(x = "k", y = NULL, color = NULL) +
    theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  ggsave(file.path(output_dir, sprintf("%s_k_diagnostics.png", m)), p, width = 8, height = 3.6, dpi = 300, bg = "white")
  
  # agreement between transforms (adjusted Rand index; 1 = identical grouping)
  prs <- combn(transforms, 2)
  ari_df <- data.frame(Month = m, n = n,
                       Pair = apply(prs, 2, paste, collapse = " vs "),
                       ARI = round(apply(prs, 2, function(p) ari(sub_df[[paste0("Cluster_", p[1])]], sub_df[[paste0("Cluster_", p[2])]])), 2))
  cat("\n  Agreement between transforms (ARI):\n"); print(ari_df[, c("Pair", "ARI")], row.names = FALSE)
  
  list(df = sub_df, summary = do.call(rbind, summ_list), ari = ari_df, k = k_used)
}

res <- lapply(months, cluster_month)
res <- res[!sapply(res, is.null)]

# -----------------------------------------------------------------------------
# 5. Save all the files 
# -----------------------------------------------------------------------------
all_df  <- do.call(rbind, lapply(res, `[[`, "df"))
keep    <- c(id_col, "Year", "Month", "Combo", "ER", "MinDist", "Conc50", paste0("Cluster_", transforms))
write.csv(all_df[, keep], file.path(output_dir, "cluster_assignment_all_transforms.csv"), row.names = FALSE)
for (tr in transforms) {
  one <- all_df[, c(id_col, "Year", "Month", "Combo")]; one$Cluster <- all_df[[paste0("Cluster_", tr)]]
  write.csv(one, file.path(output_dir, tr, sprintf("cluster_assignment_%s.csv", tr)), row.names = FALSE)
}
write.csv(do.call(rbind, lapply(res, `[[`, "summary")), file.path(output_dir, "cluster_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(res, `[[`, "ari")),     file.path(output_dir, "transform_agreement_ARI.csv"), row.names = FALSE)

