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


# scalable methods: nnSVG, HVGs

# note choice of filtering per method
spe_list <- list(
  humanDLPFC_HVGs = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_HVGs_lowFilt.rds")),
  humanDLPFC_nnSVG = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_nnSVG_lowFilt.rds")),
  humanDLPFC_SPARKX = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_SPARKX_lowFilt.rds")),
  humanDLPFC_MorI = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_MorI_lowFilt.rds"))
)


res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_nnSVG_lowFilt.rds"))),
  humanDLPFC_SPARKX = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_SPARKX_lowFilt.rds"))),
  humanDLPFC_MorI = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_MorI_lowFilt.rds"))),
  humanDLPFC_HVGs = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/spe_humanDLPFC_HVGs_lowFilt.rds")))
  )

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_MorI"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_MorI"]]), "_MorI")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_MorI$gene_id)

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



# ------------
# save results
# ------------

# save colData and spatialCoords instead of full SPE objects to save space

coldata_out <- spatialcoords_out <- list()

coldata_out$nnSVG <- colData(spe_out$spe_nnSVG)
coldata_out$SPARKX <- colData(spe_out$spe_SPARKX)
coldata_out$MorI <- colData(spe_out$spe_MorI)
coldata_out$HVGs <- colData(spe_out$spe_HVGs)

spatialcoords_out$nnSVG <- spatialCoords(spe_out$spe_nnSVG)
spatialcoords_out$SPARKX <- spatialCoords(spe_out$spe_SPARKX)
spatialcoords_out$MorI <- spatialCoords(spe_out$spe_MorI)
spatialcoords_out$HVGs <- spatialCoords(spe_out$spe_HVGs)

res_out <- list(
  coldata_out = coldata_out, 
  spatialcoords_out = spatialcoords_out
)


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/res_downstream_clustering.rds")
saveRDS(res_out, file = fn)




