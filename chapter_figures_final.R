### CHAPTER: A PRIMER ON PRE-PROCESSING, VISUALIZATION, CLUSTERING, AND PHENOTYPING OF BARCODE-BASED SPATIAL TRANSCRIPTOMICS DATA
### FROM: STATSITICAL GENOMICS
## The following code was used to produce the figures in the chapter.
## The code assumes this script is located within the following directory tree:
## 
## ./primer_spatial_transcriptomics
## │
## ├── code/
## │   ├── chapter_figures.R
## ├── figures.md
## └── data
##
## Working directory set to `./primer_spatial_transcriptomics/code`
##


## Load in Packages Needed -----------------------------------------------------
# RandomFields pulled from cran - install from github devtools::install_github("cran/RandomFields")
library('tidyverse')
library('lubridate')
library('googlesheets4')
library('ggpubr')
library('Seurat')
library('STutility')
library('spatialGE')
library('readxl')
library('preprocessCore')
library('msigdbr')
library('ComplexHeatmap')
library('BayesSpace')
#library('R.utils')
#library('hdf5r')


## Figure 1 --------------------------------------------------------------------
# Specify URL to Google Sheet from the Museum of Spatial Transcriptomics
# If not able to download, running `gs4_deauth()` may help
gs4_deauth()
mstw = 'https://docs.google.com/spreadsheets/d/1sJDb9B7AtYmfKv4-m8XR7uc3XXw_k4kGSout8cqZ8bY/edit#gid=1424019374'

# Create data frame that contains...
df_methods = read_sheet(mstw, sheet=4) %>%
  select(date_published, title, method) %>% # ROI selection
  rbind(read_sheet(mstw, sheet=5) %>%
          select(date_published, title, method)) %>% # NGS barcoding
  rbind(read_sheet(mstw, sheet=6) %>%
          select(date_published, title, method)) %>% # smFISH
  rbind(read_sheet(mstw, sheet=7) %>%
          select(date_published, title, method)) %>% # ISS
  mutate(year_pub=year(date_published)) %>%
  mutate(method=str_replace(method, " WTA| DSP", "")) %>% # Remove the first instance of " WTA" or " DSP" in the method to create "GeoMX"
  filter(method %in% c('GeoMX', 'Tomo-seq', 'Visium', 'ST', 
                       'slide-seq2', 'MERFISH', 'Molecular Cartography', 'ISS', 
                       'STARmap')) %>% # Filter down methods
  mutate(year_pub=as.factor(year_pub)) %>% # Convert year_pub to factor
  mutate(method=factor(method, levels=c("STARmap", "ISS", "ST", "Molecular Cartography", 
                                        "MERFISH", "slide-seq2", "Tomo-seq", "GeoMX", "Visium"))) %>% # Reset levels of methods
  dplyr::select(-date_published) %>% # Remove complete date (keep year of publication)
  distinct() %>% # Keep unique rows (some studies will have more than one row given diversity of tissue types)
  dplyr::count(year_pub, method) # %>% # Tally up the methods by year published
  # complete(method, year_pub) # Expand methods for those missing in years

# Creating plot of number of publications by year with stacked column plot by method
ggplot(df_methods, aes(x=year_pub, y=n, fill=method)) + 
  geom_col() +
  xlab('Year') +
  ylab('Number of publications') +
  scale_fill_muted() +
  scale_y_continuous(expand=c(0, 0)) +
  theme(panel.background=element_rect(fill="white", color="black"), 
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
# ggsave('../figures/figure_1.pdf') # Save figure


## Figure 3 --------------------------------------------------------------------
# Data from 10X Genomics (Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (Block A Section 2)
# Download and place in `data` folder located within working directory
bc_visium = Load10X_Spatial('../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/', filename='V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5')

# Using the SpatialFeaturePlot function from Seurat, plot the read counts for each spot next to H&E
ggarrange(
  Seurat::SpatialFeaturePlot(bc_visium, features="nCount_Spatial", alpha=0, image.alpha=1), # Set feature to alpha 0 and H&E to 1
  Seurat::SpatialFeaturePlot(bc_visium, features="nCount_Spatial", alpha=1, image.alpha=0, pt.size.factor=1.5), # Feature alpha to 1 and H&E alpha to 0
  ncol=2, common.legend=T, legend='bottom'
)
# ggsave('../figures/figure_3.pdf') # Save figure


## Figure 4 --------------------------------------------------------------------
bc_visium_SCT = bc_visium
bc_visium_SCT[["percent_mt"]] = PercentageFeatureSet(object=bc_visium_SCT, pattern="^MT-") # Calculate percentage of mitochondrial counts per spot
bc_visium_SCT = SCTransform(bc_visium_SCT, assay='Spatial', return.only.var.genes=F, vars.to.regress='percent_mt') # Regress out MT feature percent and normalize counts

bc_visium_log = bc_visium
bc_visium_log = NormalizeData(bc_visium_log) # Normalize the data by log transformation

# Plot raw CCNB1 (Cyclin B1) counts, log transformed, and SCT mitochondria regressed-out counts
ggarrange(
  SpatialFeaturePlot(bc_visium, features='CCNB1', image.alpha=0, pt.size.factor=1.5),
  SpatialFeaturePlot(bc_visium_log, features='CCNB1', image.alpha=0, pt.size.factor=1.5),
  SpatialFeaturePlot(bc_visium_SCT, features='CCNB1', image.alpha=0, pt.size.factor=1.5),
  ncol=3, legend='bottom', labels=c('Raw counts', 'Log-counts', 'SCTransform counts')
) 
# ggsave('../figures/figure_4.pdf') # Save figure


## Figure 5 --------------------------------------------------------------------
# Create data frame with columns containing paths to h5, position, png, and scaling factor files
files_df = data.frame(samples='../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5',
                      spotfiles='../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/spatial/tissue_positions_list.csv',
                      imgs='../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/spatial/tissue_hires_image.png',
                      json='../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/spatial/scalefactors_json.json')
bc_visium_STU = InputFromTable(infotable=files_df, platform="Visium") # Load Visium data from from data frame with file locations
bc_visium_STU = LoadImages(bc_visium_STU, time.resolve=FALSE, verbose=TRUE) # Construct sample H&E image within object
bc_visium_STU = ManualAnnotation(bc_visium_STU) # Manually annotate tissue
# Save manual annotations to text file
# write.csv(rownames_to_column(.data=bc_visium_STU@meta.data, var='barcode'), row.names=F, '../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics_STUtility_manual_annots.csv')
bc_visium_STU[["percent_mt"]] = PercentageFeatureSet(object=bc_visium_STU, pattern="^MT-") # Calculate percentage of mitochondrial counts per spot
bc_visium_STU = SCTransform(bc_visium_STU, assay='RNA', return.only.var.genes=F, vars.to.regress='percent_mt') # Regress out MT feature percent and normalize counts
colors = c('yellow', 'lightblue', 'red', 'forestgreen', 'orange', 'purple', 'pink', 'black') # Make color palette
names(colors) = unique(bc_visium_STU@meta.data$labels) # Assign cluster names to colors
ggarrange(
  ST.FeaturePlot(bc_visium_STU, features='labels', cols=colors),
  VlnPlot(bc_visium_STU, features='ESR1', group.by='labels', cols=colors),
  VlnPlot(bc_visium_STU, features='MS4A1', group.by='labels', cols=colors),
  VlnPlot(bc_visium_STU, features='CD8A', group.by='labels', cols=colors),
  ncol=2, nrow=2, common.legend=T, legend='bottom'
)
# ggsave('../figures/figure_5.pdf') # Save figure
rm(bc_visium_STU) # Clean environment

## Figure 6 --------------------------------------------------------------------
# Load data using spatialGE in the 10x Genomics breast cancer folder previously downloaded
bc_visium_spGE = STList('../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/',sample='Breast_Cancer')
bc_visium_spGE = spatialTransform(bc_visium_spGE) # Perform log transformation of the breast cancer count data
bc_visium_spGE = gene_krige(bc_visium_spGE, genes=c('ESR1', 'MS4A1')) # Perform spatial interpolation (kriring) of the expression for ESR1 and MS4A1
ggarrange(
  plotlist=append(
    plot_gene_quilt(bc_visium_spGE, genes=c('ESR1', 'MS4A1'), color_pal='BuRd'),
    plot_gene_krige(bc_visium_spGE, genes=c('ESR1', 'MS4A1'), color_pal='BuRd')),
  common.legend=T, legend='bottom', labels='AUTO'
)
# ggsave('../figures/figure_6.pdf')


## Figure 7 --------------------------------------------------------------------
# Download data from Nanostring's Spatial Organ Atlas:
# https://nanostring.com/products/geomx-digital-spatial-profiler/spatial-organ-atlas/human-kidney/
# Registration may be required before data download
# Decompress the file `hu_kidney_count_results.tar.gz` and place the resulting
# folder within then working directory under `data`
data_fp = '../data/hu_kidney_count_results/Export3_BiologicalProbeQC.xlsx'
# Read in ROI annotation and count data
roi_meta = read_excel(data_fp, sheet=1) %>% # Read Excel sheet with annotations
  janitor::clean_names() %>%
  mutate(segment_display_name=tolower(segment_display_name) %>% # Remove white spaces and other characters from ROI names
           str_replace_all("[ |]+", '_') %>%
           str_replace('panck\\+', 'pan_ck') %>%
           str_replace('cd10\\+', 'cd10')) %>% 
  mutate(morpho_marker=case_when((segment_label == 'Full ROI' | segment_label == 'Geometric Segment') ~ 'No Segmentation',
                                TRUE ~ segment_label)) # Make column with ROI segmentation info based on morphology markers
roi_counts = read_excel(data_fp, sheet=3) %>% # Read Excel sheet with count data
  janitor::clean_names() %>%
  rename_with(~gsub(pattern='^x', replacement='', x=.x)) %>% # Fix ROI names to match those in annotations
  select(c(target_name, roi_meta$segment_display_name)) %>% 
  column_to_rownames(var='target_name')

roi_quant = as.data.frame(preprocessCore::normalize.quantiles(as.matrix(roi_counts))) # Perform quantile normalization on counts
colnames(roi_quant) = colnames(roi_counts)
rownames(roi_quant) = rownames(roi_counts)

geneSets = msigdbr(species="Homo sapiens") # Get gene sets
# Pull out only biological process gene ontology
pathways_kegg = filter(geneSets, gs_subcat == "GO:BP") %>%
  split(x=.$gene_symbol, f=.$gs_name)

roi_quant_kegg = roi_quant[rownames(roi_quant) %in% pathways_kegg$GOBP_REGULATION_OF_RENAL_SYSTEM_PROCESS, ] # Extract genes related to renal regulation
roi_quant_kegg = roi_quant_kegg[, apply(roi_quant_kegg, 2, sd) != 0] # Remove genes who's standard deviation is 0
roi_quant_kegg = roi_quant_kegg[apply(roi_quant_kegg, 1, sd) > quantile(apply(roi_quant_kegg, 1, sd), 0.50), ] # Keep genes with largest 50% standard deviation
roi_meta = roi_meta[match(colnames(roi_quant_kegg), roi_meta$segment_display_name), ] # Filter ROIs in metadata that are not in count data

hm_mtx = t(scale(t(roi_quant_kegg))) # Z-scale the gene counts
roisegment = c('gray50', 'red', 'green')
names(roisegment) = as.vector(na.omit(unique(roi_meta$morpho_marker))) # Assign morphological marker names to the colors for heatmap
regioncol = as.vector(khroma::color('discreterainbow')(length(unique(roi_meta$type))))
names(regioncol) = unique(roi_meta$type) # Assign region types to the colors for heatmap
# Plot heatmap
# pdf('../figures/figure_7.pdf', width=7, height=5)
Heatmap(hm_mtx, cluster_columns=T, cluster_rows=T, show_column_names=F, show_row_dend=F,
        column_names_max_height=max_text_height(colnames(hm_mtx), gp=gpar(fontsize=500)),
        bottom_annotation=HeatmapAnnotation(df=data_frame(roi=colnames(hm_mtx), 
                                                          `ROI segment`=unlist(roi_meta$morpho_marker), 
                                                          `Spatial region`=unlist(roi_meta$type)) %>%
                                              column_to_rownames(var='roi'),
                                            col=list(`ROI segment`=roisegment, `Spatial region`=regioncol)),
        heatmap_legend_param=list(title="Scaled\nexpression")
)
# dev.off()
rm(roi_meta, roi_counts) # Clean environment

## Figure 8 --------------------------------------------------------------------
# Run the Louvain clustering workflow from Seurat
bc_visium_SCT = RunPCA(bc_visium_SCT, assay="SCT", verbose=F) %>% # Reduce dimentionality via PCA
  FindNeighbors(reduction="pca", dims=1:30) %>% # Find spot nearest neighbors from the first 30 PCs
  FindClusters(verbose=F, resolution=c(0.1, 0.8)) %>% # Run Louvain community-finding algorithm based on the nearest neighbors, using two clustering resolution levels
  RunUMAP(reduction="pca", dims=1:30) # Compute UMAP with first 30 PCs

Idents(bc_visium_SCT) = 'SCT_snn_res.0.1' # Set cluster identities for each spot

df_clusters = bc_visium_SCT@images$slice1@coordinates %>% # Extract spots' x,y locations
  rownames_to_column(var='barcode') %>% # Make row names into "barcode" column
  left_join(., bc_visium_SCT@meta.data %>%
              rownames_to_column(var='barcode'), by='barcode') %>% # Merge x,y locations with Louvain clustering assignments
  select(barcode, row, col, imagerow, imagecol, SCT_snn_res.0.1, SCT_snn_res.0.8)

bc_visium_spGE = STclust(bc_visium_spGE, ks='dtc', w=0.025) # Apply STclust algorithm with spatial weight of 0.025 and dynamic tree cut to split hierarchical clusteing dendrogram
df_clusters = df_clusters %>% # Merge Louvain results with STclust results
  left_join(., bc_visium_spGE@st_clusters$Breast_Cancer_spw0.025, by=c('row'='ypos', 'col'='xpos')) %>%
  dplyr::rename(libname_dtc=libname, WCluster_dtc=WCluster)

bc_visium_spGE = STclust(bc_visium_spGE, ks=15, w=0.025) # Apply STclust algorithm with spatial weight of 0.025 and split hierarchical clusteing dendrogram into 15 clusters (k=15)
df_clusters = df_clusters %>% # Merge STclust results with previous clustering outputs
  left_join(., bc_visium_spGE@st_clusters$Breast_Cancer_spw0.025, by=c('row'='ypos', 'col'='xpos')) %>%
  dplyr::rename(libname_k15=libname, WCluster_k15=WCluster)

bc_visium_BSP = readVisium('../data/Breast_Cancer_Ductal_Carcinoma_10XGenomics/') # Load in breast cancer Visium data to be used in BayesSpace
bc_visium_BSP = spatialPreprocess(bc_visium_BSP, platform="Visium", n.PCs=10, n.HVGs=2000) # Process and transform data prior to clustering (Use first 10 PCs and 2000 most variable genes)
set.seed(149)
bc_visium_BSP = spatialCluster(bc_visium_BSP, q=15, platform="Visium", d=10, # Run BayesSpace clustering
                               init.method="mclust", model="t", gamma=2,
                               nrep=1000, burn.in=100,
                               save.chain=TRUE)

# Merge previous clustering results with clusters identified using BayesSpace
df_clusters = df_clusters %>%
  left_join(., colData(bc_visium_BSP) %>%
              as.data.frame() %>%
              select(spot, spatial.cluster), by=c('barcode'='spot')) %>%
  dplyr::rename(bayessp=spatial.cluster) %>%
  mutate(bayessp=as.factor(.$bayessp))

# Run BayesSpace sub-spot resolution clustering algorithm
bc_visium_BSP_enhanced = spatialEnhance(bc_visium_BSP, q=15, platform="Visium", d=10,
                                        model="t", gamma=2,
                                        jitter_prior=0.3, jitter_scale=3.5,
                                        nrep=1000, burn.in=100,
                                        save.chain=TRUE)

# Pull out sub-spot clustering assignments
df_clusters_enh = colData(bc_visium_BSP_enhanced) %>%
  as.data.frame() %>%
  select(imagerow, imagecol, spatial.cluster) %>%
  mutate(bayessp=as.factor(.$spatial.cluster))

# Create list of plots with clustering solutions to arrange together
cluster_p = list()
cluster_p[['seurat_0.1']] = ggplot(df_clusters, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=SCT_snn_res.0.1), size=0.4) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()
cluster_p[['seurat_0.8']] = ggplot(df_clusters, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=SCT_snn_res.0.8), size=0.4) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()
cluster_p[['STclustDTC']] = ggplot(df_clusters, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=WCluster_dtc), size=0.4) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()
cluster_p[['STclustk15']] = ggplot(df_clusters, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=WCluster_k15), size=0.4) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()
cluster_p[['bayesspace']] = ggplot(df_clusters, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=bayessp), size=0.4) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()
cluster_p[['bayesspaceEnh']] = ggplot(df_clusters_enh, aes(x=imagecol, y=imagerow)) +
  geom_point(aes(color=bayessp), size=0.01) +
  scale_y_reverse() +
  khroma::scale_color_discreterainbow() +
  theme_void() +
  theme(legend.position="none") +
  coord_equal()

# NOTE: Panels G and H of Figure 8 were generated using the Python code in "chapter_figures_spaGCN.py"

ggarrange(plotlist=cluster_p, ncol=2, nrow=4,  align='hv', labels='AUTO')
# ggsave('../figures/figure_8.pdf', width=5)
rm(bc_visium_spGE, bc_visium_BSP, bc_visium_BSP_enhanced) # Clean environment

## Figure 9 --------------------------------------------------------------------
# Code based on SPOTlight tutorial developed by the algorithm's authors
# NOTE: Must close and re-open R because of conflicts with loaded packages and `scran::scoreMarkers`
library('TabulaMurisSenisData')
library('TENxVisiumData')
library('SPOTlight')
sce = TabulaMurisSenisDroplet(tissues="Kidney")$Kidney # Download single cell data set for mouse kidney (reference data set)
spe = MouseKidneyCoronal() # Download Visium spatial data for mouse kidney (data to be deconvolved)
rownames(spe) = rowData(spe)$symbol

sce = subset(sce, , age == "18m") # Subset single cell data to only cells from the 18-month age samples
sce = subset(sce, , !free_annotation %in% c("nan", "CD45")) # Remove cells with CD45 and missing annotation
sce = scuttle::logNormCounts(sce) # Log transform the counts
dec = scran::modelGeneVar(sce) # Model mean-variance relationship for each gene
hvg = scran::getTopHVGs(dec, n=3000) # Select the top 3000 most variable genes

colLabels(sce) = colData(sce)$free_annotation # Set the cell type annotation
genes = !grepl(pattern="^Rp[l|s]|Mt", x=rownames(sce)) # Remove ribosomal and mitochondrial genes from single-cell reference data
mgs = scran::scoreMarkers(sce, subset.row=genes) # Score genes to use in cell type assignment

mgs_fil = lapply(names(mgs), function(i) { #for each of the cell types
  x = mgs[[i]]
  # Filter and keep relevant marker genes, those with AUC > 0.8
  x = x[x$mean.AUC > 0.8, ]
  # Sort the genes from highest to lowest weight
  x = x[order(x$mean.AUC, decreasing = TRUE), ]
  # Add gene and cluster id to the data frame
  x$gene = rownames(x)
  x$cluster = i
  data.frame(x)
})
mgs_df = do.call(rbind, mgs_fil) # Paste together all data frames in the list 

# Create a list of cell types that contain the index of each
idx = split(seq(ncol(sce)), sce$free_annotation) 
# Downsample to at most 20 cells per cell type and subset (to speed up SPOTlight deconvolution)
n_cells = 100
cs_keep = lapply(idx, function(i) {
  n = length(i)
  if (n < n_cells)
    n_cells = n
  sample(i, n_cells)
})
sce = sce[, unlist(cs_keep)] # Select columns in SingleCellExperiment object corresponding to the subsampled cells

# Run SPOTlight deconvolution
res = SPOTlight(
  x=sce,
  y=spe,
  groups=sce$free_annotation,
  mgs=mgs_df,
  hvg=hvg,
  weight_id="mean.AUC",
  group_id="cluster",
  gene_id="gene")

mat = res$mat # Extract deconvolved scores for each spot

ct = colnames(mat) # Get cell types in data set
mat[mat < 0.1] = 0 # Assign zero to scores less than 0.1 (less cell types to plot in scatterpie)
pal = khroma::color('discreterainbow')(length(ct)) # Get color palette
names(pal) = ct # Name the cell type vector with the colors

library('ggplot2')
plotSpatialScatterpie( # Plot the deconvolution results for each spot
  x=spe,
  y=mat,
  cell_types=ct,
  img=F,
  scatterpie_alpha=1,
  pie_scale=0.4) +
  scale_fill_manual(values=pal, breaks=names(pal))

# ggsave('../figures/figure_9.pdf')

