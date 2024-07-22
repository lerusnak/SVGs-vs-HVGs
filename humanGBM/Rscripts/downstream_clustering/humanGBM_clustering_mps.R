######################
#     human GBM      #
#   dim reduction    #
#   & clustering     #
#   Meta Programs    #
######################


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

# GBM meta data - for cluster labels

visium_metadata <- read.csv("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/visium_metadata.csv")

spot_labels <- visium_metadata %>%
  select(spot_id, sample, mp, layer) %>%
  filter(sample == "MGH258") 
dim(spot_labels)


# spots included in GBM clustering analysis
spots <- spot_labels$spot_id
length(spots)


# scalable methods: nnSVG, HVGs

# note choice of filtering per method

humanGBM_HVGs <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_HVGs.rds"))
humanGBM_nnSVG <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_nnSVG.rds"))
humanGBM_SPARKX <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_SPARKX.rds"))
humanGBM_MorI <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_MorI.rds"))

humanGBM_HVGs <- humanGBM_HVGs[,colnames(humanGBM_HVGs) %in% spots]
dim(humanGBM_HVGs)

humanGBM_nnSVG <- humanGBM_nnSVG[,colnames(humanGBM_nnSVG) %in% spots]
dim(humanGBM_nnSVG)

humanGBM_SPARKX<- humanGBM_SPARKX[,colnames(humanGBM_SPARKX) %in% spots]
dim(humanGBM_SPARKX)

humanGBM_MorI <- humanGBM_MorI[,colnames(humanGBM_MorI) %in% spots]
dim(humanGBM_MorI)

spe_list <- list(humanGBM_HVGs = humanGBM_HVGs, 
                 humanGBM_nnSVG = humanGBM_nnSVG, 
                 humanGBM_SPARKX = humanGBM_SPARKX, 
                 humanGBM_MorI = humanGBM_MorI)

res_list <- list(humanGBM_HVGs = rowData(humanGBM_HVGs), 
                 humanGBM_nnSVG = rowData(humanGBM_nnSVG),
                 humanGBM_SPARKX = rowData(humanGBM_SPARKX),
                 humanGBM_MorI = rowData(humanGBM_MorI))


# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanGBM_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanGBM_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanGBM_MorI"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_MorI"]]), "_MorI")[-(1:2)]
colnames(res_list[["humanGBM_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$humanGBM_HVGs$gene_id %in% res_list$humanGBM_nnSVG$gene_id)
table(res_list$humanGBM_HVGs$symbol %in% res_list$humanGBM_SPARKX$symbol)
table(res_list$humanGBM_HVGs$symbol %in% res_list$humanGBM_MorI$symbol)

spe_out <- list()



# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanGBM_HVGs$rank_HVGs <= 1000)
top <- res_list$humanGBM_HVGs$gene_id[ix]

spe <- spe_list$humanGBM_HVGs

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(1234)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(1234)
g <- buildSNNGraph(spe, k = 8, use.dimred = "PCA")
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
ix <- which(res_list$humanGBM_nnSVG$rank_nnSVG <= 1000)
top <- res_list$humanGBM_nnSVG$gene_id[ix]

spe <- spe_list$humanGBM_nnSVG

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(1234)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(1234)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(1234)
g <- buildSNNGraph(spe, k = 7, use.dimred = "PCA")
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
ix <- which(res_list$humanGBM_SPARKX$rank_SPARKX <= 1000)
genes_subset <- res_list$humanGBM_SPARKX[ix,]
top <- rownames(genes_subset)

spe <- spe_list$humanGBM_SPARKX

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(1234)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(1234)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(1234)
g <- buildSNNGraph(spe, k = 8, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_SPARKX <- spe



# --------------------------------
# downstream clustering: Moran's I
# --------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# Mor I: top 1000 SVGs
ix <- which(res_list$humanGBM_MorI$rank_MorI <= 1000)
genes_subset <- res_list$humanGBM_MorI[ix,]
top <- rownames(genes_subset)

spe <- spe_list$humanGBM_MorI

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(1234)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(1234)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(1234)
g <- buildSNNGraph(spe, k = 8, use.dimred = "PCA")
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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/res_downstream_clustering_mps.rds")
saveRDS(res_out, file = fn)




