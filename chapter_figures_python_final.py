### CHAPTER: A PRIMER ON PRE-PROCESSING, VISUALIZATION, CLUSTERING, AND PHENOTYPING OF BARCODE-BASED SPATIAL TRANSCRIPTOMICS DATA
### FROM: STATSITICAL GENOMICS
## The following Python code was used to produce panels G and H of Figure 8 in the chapter.
## It also assumes that the data used in the R script `chapter_figures_final.R` has been placed in the `data` subfolder.
## The code assumes this script is located within the following directory tree:
##
## ./primer_spatial_transcriptomics
## │
## ├── code/
## │   ├── chapter_figures_final.R
## │   └── chapter_figures_python_final.py
## ├── figures.md
## └── data
##
## Working directory set to `./primer_spatial_transcriptomics/code`
##


## Import modules ------------------------------------------------------------------------------------------------------
from pathlib import Path
from scanpy import read_10x_h5
import re
import cv2
import pandas as pd
import SpaGCN as spg
import scanpy as sc

# Get working directory
wd = re.sub('code', '', str(Path.cwd()))


## Figure 8 (Panels G and H) -------------------------------------------------------------------------------------------
# Get paths to HDF5, spatial coordinates, and image files
# NOTE: Files must be placed in the `data` subfolder
# NOTE: Code based on SpaGCN's tutorial created by the package's authors
expr_path = wd+"data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/V1_Breast_Cancer_Block_A_Section_2_raw_feature_bc_matrix.h5"
spat_path = wd+"data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/spatial/tissue_positions_list.csv"
image_path = wd+"data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/spatial/tissue_hires_image.png"

# Read in expression data
bc_visium_py = read_10x_h5(expr_path)
# Read in spatial coordinates
spatial_df = pd.read_csv(spat_path, sep=",", header=None, na_filter=False, index_col=0)
# Save spatial data in `anndata` object
# Following the tutorial, x represents rows of the 10X Visium array in the tissue_positions_list.csv file
bc_visium_py.obs["in_tissue"] = spatial_df[1]
bc_visium_py.obs["x_array"] = spatial_df[2]
bc_visium_py.obs["y_array"] = spatial_df[3]
bc_visium_py.obs["x_pixel"] = spatial_df[4]
bc_visium_py.obs["y_pixel"] = spatial_df[5]
# Scale pixel coordinates using factor (0.08250825) in json file produced by the Visium platform
bc_visium_py.obs["x_pixel_sc"] = bc_visium_py.obs["x_pixel"] * 0.08250825
bc_visium_py.obs["y_pixel_sc"] = bc_visium_py.obs["y_pixel"] * 0.08250825
# Select only spots covered with tissue (binary 0:No tissue over spot ;  1:Tissue over spot)
bc_visium_py = bc_visium_py[bc_visium_py.obs["in_tissue"] == 1]
# Force gene names to strings
bc_visium_py.var['genename'] = bc_visium_py.var.index.astype('str')

# Read image file
img = cv2.imread(image_path)
# Get coordinate data from anndata
x_array = bc_visium_py.obs['x_array'].tolist()
y_array = bc_visium_py.obs['y_array'].tolist()
# Convert coordinates to integers given that scaling create float numbers
x_pixel = bc_visium_py.obs['x_pixel_sc'].tolist()
x_pixel = [int(i) for i in x_pixel]
y_pixel = bc_visium_py.obs['y_pixel_sc'].tolist()
y_pixel = [int(i) for i in y_pixel]

# Make sure that spots are aligned to image
# This part of the code generates an image of the H&E stained tissue, overlaid with squares re[resenting the Visium array spots
# The squares must match the general shape of the tissue slice
# The `res` variable controls the size of the square
res = 5
img_new = img.copy()
for i in range(len(x_pixel)):
    x = x_pixel[i]
    y = y_pixel[i]
    img_new[int(x-res):int(x+res), int(y-res):int(y+res), : ] = 255
# Image is written to the `code` sub folder for review
cv2.imwrite('./check_spot_alignment.jpg', img_new)

# Run SpaGCN WITHOUT histology weight (Panel G)
# This step predicts the clusters and then stores the cluster assignments within anndata
bc_visium_py.obs["pred_s0"] = spg.detect_spatial_domains_ez_mode(bc_visium_py, img, x_array, y_array, x_pixel, y_pixel, n_clusters=15, histology=True, s=0, b=49, p=0.5, r_seed=100, t_seed=100, n_seed=100)
bc_visium_py.obs["pred_s0"] = bc_visium_py.obs["pred_s0"].astype('category')

# Plot clusters WITHOUT histology weight (Panel G)
# Specify colors for plotting
colors = ["#D1BBD7", "#AE76A3", "#882E72", "#1965B0", "#5289C7", "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7F056", "#F6C141", "#F1932D", "#E8601C", "#DC050C", "#72190E"]
# spg.plot_spatial_domains_ez_mode(bc_visium_py, domain_name="pred_s0", x_name="y_pixel", y_name="x_pixel", plot_color=colors, size=150000/bc_visium_py.shape[0], show=False, save_dir="../figures/figure_8_panel_G.pdf")

# Run SpaGCN WITH histology weight (Panel H)
# This step predicts the clusters and then stores the cluster assignments within anndata
bc_visium_py.obs["pred_s2"] = spg.detect_spatial_domains_ez_mode(bc_visium_py, img, x_array, y_array, x_pixel, y_pixel, n_clusters=15, histology=True, s=2, b=49, p=0.5, r_seed=100, t_seed=100, n_seed=100)
bc_visium_py.obs["pred_s2"] = bc_visium_py.obs["pred_s2"].astype('category')

# Plot clusters WITH histology weight (Panel H)
# spg.plot_spatial_domains_ez_mode(bc_visium_py, domain_name="pred_s2", x_name="y_pixel", y_name="x_pixel", plot_color=colors, size=150000/bc_visium_py.shape[0], show=False, save_dir="../figures/figure_8_panel_H.pdf")

