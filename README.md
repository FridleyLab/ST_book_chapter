# A primer on pre-processing, visualization, clustering, and phenotyping of barcode-based spatial transcriptomics
## By: Oscar Ospina, Alex Soupir, and Brooke L. Fridley <br/> From: Statistical Genomics

This repository contains the code to create the figures in the book chapter "A primer on 
pre-processing, visualization, clustering, and phenotyping of barcode-based spatial transcriptomics",
part of the book (under preparation) entitled "Statistical Genomics" (Eds: Drs. Brooke L. Fridley and
Xuefeng Wang).

Data used to create the figures in this chapter are not included in this repository and will
need to be downloaded by the reader from the following sources:

### Figure 1:
Museum of Spatial Transcriptomics (by Moses and Pachter, 2022): Excel spreadsheet containing 
metadata on published articles and preprints. The code contains commands to directly download 
and read the data within R via `googlesheets4`. [LINK TO PAPER DESCRIBING THE DATABASE](https://www.nature.com/articles/s41592-022-01409-2).

### Figures 3, 4, 5, 6, and 8:
Ductal breast carcinoma 10X Visium sample, available at the 10X data repository. Before
download, registration to 10X Genomics may be required. To download the data, The reader can 
query "Human Breast Cancer Block A Section 2", select "Spatial Gene Expression", and
then, select "Space Ranger 1.1.0". Clicking on the query result will take the reader to the sample
download webpage. The code in this chapter uses the barcode matrices and the spatial imaging
data. [LINK TO 10X GENOMICS DATA REPOSITORY](https://www.10xgenomics.com/resources/datasets).

### Figure 7:
Human kidney GeoMx sample, available at Nanostring's Spatial Organ Atlas website. Before
download, registration to Nanostring may be required. To download the data, the reader can 
click on "Data Downloads" and then click on "Count Results". This link will produce a compressed
`.gz` file, which must be decompressed by the reader. The resulting folder contains the Excel
spreadsheet with ROI metadata and gene counts necessary to generate the chapter's figure.
[LINK TO NANOSTRING SPATIAL ORGAN ATLAS](https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/human-kidney/).

### Figure 9:
Mouse kidney single-cell and spatial data sets from the R libraries `TabulaMurisSenisData`
(Soneson et al. 2022) [GITHUB](https://github.com/fmicompbio/TabulaMurisSenisData) and 
`TENxVisiumData` (Crowell, 2022) [GITHUB](https://github.com/helenalc/TENxVisiumData). 
The code contains commands to directly download and read the data sets within R.

