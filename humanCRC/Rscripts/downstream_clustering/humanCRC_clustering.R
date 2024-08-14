######################
#     human CRC      #
#   dim reduction    #
#   & clustering     #
#      LAYERS        #
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

# scalable methods: SPARK-X, SpatialDE2, HVGs

### Load spatial experiment objects for all VG methods ###

humanCRC_HVGs <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_HVGs.rds"))
humanCRC_SPARKX <- readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_SPARKX.rds"))


spe_list <- list(humanCRC_HVGs = humanCRC_HVGs, 
                 humanCRC_SPARKX = humanCRC_SPARKX)

res_list <- list(humanCRC_HVGs = rowData(humanCRC_HVGs), 
                 humanCRC_SPARKX = rowData(humanCRC_SPARKX))


# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanCRC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanCRC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanCRC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanCRC_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$humanCRC_HVGs$symbol %in% res_list$humanCRC_SPARKX$symbol)

spe_out <- list()


# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanCRC_HVGs$rank_HVGs <= 1000)
genes_subset <- res_list$humanCRC_HVGs[ix,]
top <- rownames(genes_subset)

spe <- spe_list$humanCRC_HVGs

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



# ------------------------------
# downstream clustering: SPARK-X
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# SPARK-X: top 1000 SVGs
ix <- which(res_list$humanCRC_SPARKX$rank_SPARKX <= 1000)
genes_subset <- res_list$humanCRC_SPARKX[ix,]
top <- rownames(genes_subset)

spe <- spe_list$humanCRC_SPARKX

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
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_SPARKX <- spe



# ------------
# save results
# ------------

# save colData and spatialCoords instead of full SPE objects to save space

coldata_out <- spatialcoords_out <- list()

coldata_out$SPARKX <- colData(spe_out$spe_SPARKX)
coldata_out$HVGs <- colData(spe_out$spe_HVGs)

spatialcoords_out$SPARKX <- spatialCoords(spe_out$spe_SPARKX)
spatialcoords_out$HVGs <- spatialCoords(spe_out$spe_HVGs)


res_out <- list(
  coldata_out = coldata_out, 
  spatialcoords_out = spatialcoords_out
)


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/res_downstream_clustering_layers.rds")
saveRDS(res_out, file = fn)

