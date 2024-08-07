############################
#        human DLPFC       #
# dimensionality reduction #
#        & clustering      #
############################


#############
# packages  #
#############

library(here)
library(SpatialExperiment)
library(tidyverse)
library(scater)
library(scran)

###############
#  load data  #
###############

### Create spe for SpatialDE2 data ###

# Spatial experiment object to use for SpatialDE2 output
sde2_spe <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_HVGs_lowFilt.rds"))

# Load run SpatialDE2 python output
sde2.py.out <- read.csv(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/humanDLPFC_SpatialDE2.csv"))

# Coerce gene_id to rownames in SpatialDE2 python output (rowData)
sde2.rowdat <- sde2.py.out %>% column_to_rownames(var  ="X")

# Sort SDE2 rowdata to match spe 
sde2.rowdat.sorted <- sde2.rowdat[rownames(sde2_spe), ]

# check that spe and SDE2.rowdat.sorted have same number of features,
# and all rownames of features match
dim(sde2.rowdat.sorted)
dim(rowData(sde2_spe))
all(rownames(rowData(sde2_spe)) == sde2.rowdat.sorted$gene_id)

# Overwrite spe rowdata with sorted SpatialDE2 python output
rowData(sde2_spe) <- sde2.rowdat.sorted


### load in spatial experiments from all VG methods ###

spe_list <- list(
  humanDLPFC_HVGs = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_HVGs_lowFilt.rds")),
  humanDLPFC_nnSVG = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_nnSVG_lowFilt.rds")),
  humanDLPFC_SPARKX = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_SPARKX_lowFilt.rds")),
  humanDLPFC_MorI = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_MorI_lowFilt.rds")),
  humanDLPFC_SpatialDE2 = sde2_spe
)


res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_nnSVG_lowFilt.rds"))),
  humanDLPFC_SPARKX = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_SPARKX_lowFilt.rds"))),
  humanDLPFC_MorI = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_MorI_lowFilt.rds"))),
  humanDLPFC_HVGs = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_HVGs_lowFilt.rds"))),
  humanDLPFC_SpatialDE2 = rowData(sde2_spe)
  )

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_MorI"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_MorI"]]), "_MorI")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["humanDLPFC_SpatialDE2"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SpatialDE2"]]), "_SpatialDE2")[-(1:2)]


# note filtering per method

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_MorI$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_SpatialDE2$gene_id)

spe_out <- list()



# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanDLPFC_HVGs$rank_HVGs <= 1000)
top <- res_list$humanDLPFC_HVGs$gene_id[ix]

spe <- spe_list$humanDLPFC_HVGs

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 25, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_HVGs <- spe





# ----------------------------
# downstream clustering: nnSVG
# ----------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# nnSVG: top 1000 SVGs
ix <- which(res_list$humanDLPFC_nnSVG$rank_nnSVG <= 1000)
top <- res_list$humanDLPFC_nnSVG$gene_id[ix]

spe <- spe_list$humanDLPFC_nnSVG

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_nnSVG <- spe



# ------------------------------
# downstream clustering: SPARK-X
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# SPARK-X: top 1000 SVGs
ix <- which(res_list$humanDLPFC_SPARKX$rank_SPARKX <= 1000)
top <- res_list$humanDLPFC_SPARKX$gene_id[ix]

spe <- spe_list$humanDLPFC_SPARKX

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 15, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_SPARKX <- spe




# ------------------------------
# downstream clustering: MoransI
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# MorI: top 1000 SVGs
ix <- which(res_list$humanDLPFC_MorI$rank_MorI <= 1000)
top <- res_list$humanDLPFC_MorI$gene_id[ix]

spe <- spe_list$humanDLPFC_MorI

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 6, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_MorI <- spe




# ----------------------------------
# downstream clustering: SpatialDE2
# ----------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 SpatialDE2 VGs
ix <- which(res_list$humanDLPFC_SpatialDE2$rank_SpatialDE2 <= 1000)
top <- res_list$humanDLPFC_SpatialDE2$gene_id[ix]

spe <- spe_list$humanDLPFC_SpatialDE2

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 6, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_SpatialDE2 <- spe



# ------------
# save results
# ------------

# save colData and spatialCoords instead of full SPE objects to save space

coldata_out <- spatialcoords_out <- list()

coldata_out$nnSVG <- colData(spe_out$spe_nnSVG)
coldata_out$SPARKX <- colData(spe_out$spe_SPARKX)
coldata_out$MorI <- colData(spe_out$spe_MorI)
coldata_out$SpatialDE2 <- colData(spe_out$spe_SpatialDE2)
coldata_out$HVGs <- colData(spe_out$spe_HVGs)

spatialcoords_out$nnSVG <- spatialCoords(spe_out$spe_nnSVG)
spatialcoords_out$SPARKX <- spatialCoords(spe_out$spe_SPARKX)
spatialcoords_out$MorI <- spatialCoords(spe_out$spe_MorI)
spatialcoords_out$SpatialDE2 <- spatialCoords(spe_out$spe_SpatialDE2)
spatialcoords_out$HVGs <- spatialCoords(spe_out$spe_HVGs)

res_out <- list(
  coldata_out = coldata_out, 
  spatialcoords_out = spatialcoords_out
)


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/res_downstream_clustering.rds")
saveRDS(res_out, file = fn)




