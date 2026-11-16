# 1. Overview

This repository contains the scripts and data used to (i) build an
updated bioavailable strontium isoscape for contiguous North America
(plant samples only), (ii) align the hydrogen, sulfur, and strontium
isoscapes, (iii) run geographic assignments for each sample, (iv)
calculate enrichment ratios and classify samples as putative migrants or
locals, (v) cluster putative migrants by month, (vi) simulate wind
trajectories, and (vii) produce the figures. All analyses are in R. The
seven numbered folders follow the order of the workflow. Outputs of one
folder are inputs to the next, so run the folders in numerical order.

# 2. Data availability

## Included in this repository

All data needed to run the scripts, except for the items listed below,
are included.

## Not shared

- Part of the strontium data used for the updated isoscape is pending
  publication and cannot be shared yet. It will be made public once
  published. Without it, Folder 1 cannot reproduce the published
  isoscape exactly. The resulting isoscape files are provided in Folder
  1 and 3 so that Folders 2 to 7 can still be run.

## Obtain separately

- Rasters used to build the strontium isoscape: available through the
  references listed in the supplementary material.

- GDAS1 meteorological files (about 42.3 GB), used in Folder 6:
  available from the NOAA website
  (https://www.ready.noaa.gov/data/archives/gdas1/). Place them in the
  specified folder before running the Folder 6 scripts.

# 3. Folder contents

## 1_Strontium_NA_isoscape

Builds the updated bioavailable strontium isoscape using Random Forest.
The output is a TIF file with the mean predicted isoscape and its
uncertainty across contiguous North America. The strontium data are
partly available from the references cited in the script; the remainder
is pending publication (see Section 2).

| **Script / file**    | **Input**                      | **Output**                           |
|----------------------|--------------------------------|--------------------------------------|
| Strontium_isoscape.R | Strontium data and rasters     | Sr.tif (mean and uncertainty layers) |
| extract_function.R   | Called by Strontium_isoscape.R | Helper function from C. Bataille.    |

## 2_Isoscapes_alignment

Obtains the hydrogen isoscape from the data included in the assignR
package and aligns it to the sulfur isoscape from Dargent et al., 2026
and to the strontium isoscape from Folder 1. "Aligned" here means
resampled to a common resolution, extent, and projection. The aligned
isoscapes are saved in Folder 3.

| **Script / file**  | **Input**                                                                                             | **Output**                           |
|--------------------|-------------------------------------------------------------------------------------------------------|--------------------------------------|
| Isoscapes_aligment | assignR hydrogen isoscape; sulfur isoscape from Dargent et al. 2026; strontium isoscape from Folder 1 | Three aligned isoscapes (H, S, Sr) i |

## 3_Assignments

Contains one R script, one Excel data file with the isotope values of
each sample, and the three aligned isoscapes (hydrogen, sulfur,
strontium). The script uses the assignR package to run geographic
assignments for every sample, once for each combination of isotopes. It
creates one output folder per combination. For each sample, each folder
contains four files:

- continuous posterior probability surface (TIFF);

- the same surface normalized (PNG);

- binary surface (TIFF), marking the area with at least 2:1 odds ratio,
  i.e., the cells that make up the highest-probability region where the
  odds of origin are 2:1 or greater;

- the binary surface as a PNG.

The continuous and binary TIFF files are used in Folders 4 to 7. The
script also combines the assignR output on weighted distance and bearing
from the sampling site to the probability-weighted origin into one
table.

<table>
<colgroup>
<col style="width: 31%" />
<col style="width: 34%" />
<col style="width: 34%" />
</colgroup>
<thead>
<tr class="header">
<th><strong>Script / file</strong></th>
<th><strong>Input</strong></th>
<th><strong>Output</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1_final_assignment_aligned</td>
<td><p>data.xlsx;</p>
<p>three isoscapes (.tif)</p></td>
<td>One folder per isotope combination with four files per sample;
distance and bearing table file</td>
</tr>
</tbody>
</table>

## 4_Enrichment_ratio

Calculates, for each sample, how much of the continuous posterior
probability falls inside the sampling zones, relative to the area those
zones cover. The enrichment ratio is the proportion of probability mass
inside the sampling zones ÷ the proportion of the total area covered by
the sampling zones. A ratio of 1 or higher means the sample is
classified as having a putative origin inside the sampling zones. The
script reads the sampling-zone point file and the folder of continuous
posterior probability TIFF files and writes an Excel file with each
sample's classification.

| **Script / file**  | **Input**                                                        | **Output**                                                                              |
|--------------------|------------------------------------------------------------------|-----------------------------------------------------------------------------------------|
| 1_Enrichment_ratio | point.xlsx; continuous posterior probability TIFFs from Folder 3 | local_migratory_classification.xlsx with enrichment ratio and classification per sample |

## 5_Cluster

Takes the classification from Folder 4 and, using the same continuous
posterior probability TIFFs, groups the putative migrants by the month
recorded in the enrichment results. It also compares clustering results
with and without transformation to show how the transformation affects
the clusters.

| **Script / file** | **Input**                                                                       | **Output**                            |
|-------------------|---------------------------------------------------------------------------------|---------------------------------------|
| 1_cluster_month   | Classification Excel file from Folder 4; continuous posterior probability TIFFs | cluster_assignment_all_transforms.csv |

## 6_Trajectory

Contains three scripts that must be run in order, plus one RDS file
created by the second script. GDAS1 meteorological files are required
(see Section 2).

| **Script / file**                    | **Input**                                                 | **Output**                                                                                    |
|--------------------------------------|-----------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| 1_simulated_trajectories             | GDAS1 files; sample information                           | Simulated trajectory for each sample, using the settings described in the manuscript.         |
| 2_Filter_warm_trajectories           | Trajectories from script 1                                | Warm trajectories only (RDS file: retained_trajectories).                                     |
| 3_Plot_warm_trajectories_per_cluster | RDS file from script 2; cluster assignments from Folder 5 | Warm trajectories grouped by the clusters from Folder 5: <span class="mark">\[output\]</span> |

## 7_Figures

Each script plots putative migrants and locals together on the same map,
using the binary TIFF files for each group within a month. The panels
created here were later assembled into the final figures in CorelDRAW.

<table>
<colgroup>
<col style="width: 26%" />
<col style="width: 36%" />
<col style="width: 36%" />
</colgroup>
<thead>
<tr class="header">
<th><strong>Script / file</strong></th>
<th><strong>Input</strong></th>
<th><strong>Output</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><p>1_Migratory</p>
<p>2_Locals</p></td>
<td>Binary TIFFs from Folder 3; classification from Folder 4</td>
<td>Figure panels, file format PNG and SVG.</td>
</tr>
</tbody>
</table>

# Glossary

| **Term**                          | **Meaning**                                                                                                                            |
|-----------------------------------|----------------------------------------------------------------------------------------------------------------------------------------|
| Isoscape                          | A map of the predicted spatial variation of an isotope ratio.                                                                          |
| Posterior probability surface     | For one sample, the probability that each map cell is the origin, given the sample's isotope values (assignR output).                  |
| Putative migrant / putative local | Classification of a sample based on its enrichment ratio (Folder 4). "Putative" because the origin is inferred, not directly observed. |
