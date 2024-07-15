############################
#       mouse Embryo       #
# dimensionality reduction #
#       & clustering       #
############################


#############
# packages  #
#############

library(here)
library(SpatialExperiment)
library(tidyverse)
library(scater)
library(scran)
library(readxl)

###############
#  load data  #
###############

mouseEmbryo_HVGs_spe <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_HVGs.rds"))
dim(mouseEmbryo_HVGs_spe)
mouseEmbryo_nnSVG_spe <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_nnSVG.rds"))
dim(mouseEmbryo_nnSVG_spe)

## Top ranked genes for each tissue in mouse embryo slide seq paper

genes_file <- "/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/data/41588_2023_1435_MOESM4_ESM.xlsx"
genes_table <- read_excel(path = genes_file, sheet = "Supplementary Table 4", skip = 3)

genes_E8.5 <- genes_table[-1, grepl("E8.5", names(genes_table))]
genes_E8.5 <- genes_E8.5[-c(21:38),]

genes_E8.5 <- unlist(genes_E8.5)


# Subset with top ranked genes (based on paper that used SPARK)

HVGs_spe_sub <- mouseEmbryo_HVGs_spe[rownames(mouseEmbryo_HVGs_spe) %in% genes_E8.5]
dim(HVGs_spe_sub)

nnSVG_spe_sub <- mouseEmbryo_nnSVG_spe[rownames(mouseEmbryo_nnSVG_spe) %in% genes_E8.5]
dim(nnSVG_spe_sub)

# scalable methods: nnSVG, HVGs

# note choice of filtering per method
spe_list <- list(
  mouseEmbryo_HVGs = HVGs_spe_sub,
  mouseEmbryo_nnSVG = nnSVG_spe_sub
)

res_list <- list(
  mouseEmbryo_nnSVG = rowData(HVGs_spe_sub),
  mouseEmbryo_HVGs = rowData(nnSVG_spe_sub))

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseEmbryo_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["mouseEmbryo_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)

spe_out <- list()



# ---------------------------
# downstream clustering: HVGs
# ---------------------------

spe <- spe_list$mouseEmbryo_HVGs

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(2)
spe <- runPCA(spe)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 70, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_HVGs <- spe





# ----------------------------
# downstream clustering: nnSVG
# ----------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

spe <- spe_list$mouseEmbryo_nnSVG

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123)
spe <- runPCA(spe)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 70, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_nnSVG <- spe



# ------------
# save results
# ------------

# save colData and spatialCoords instead of full SPE objects to save space

coldata_out <- spatialcoords_out <- list()

coldata_out$nnSVG <- colData(spe_out$spe_nnSVG)
coldata_out$HVGs <- colData(spe_out$spe_HVGs)

spatialcoords_out$nnSVG <- spatialCoords(spe_out$spe_nnSVG)
spatialcoords_out$HVGs <- spatialCoords(spe_out$spe_HVGs)

res_out <- list(
  coldata_out = coldata_out, 
  spatialcoords_out = spatialcoords_out
)


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/res_downstream_clustering_subset.rds")
saveRDS(res_out, file = fn)




