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

###############
#  load data  #
###############


# scalable methods: nnSVG, HVGs

# note choice of filtering per method
spe_list <- list(
  mouseEmbryo_HVGs = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_HVGs.rds")),
  mouseEmbryo_nnSVG = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_nnSVG.rds"))
)

res_list <- list(
  mouseEmbryo_nnSVG = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_nnSVG.rds"))),
  mouseEmbryo_HVGs = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_HVGs.rds"))))

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseEmbryo_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["mouseEmbryo_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)

spe_out <- list()



# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$mouseEmbryo_HVGs$rank_HVGs <= 1000)
top <- res_list$mouseEmbryo_HVGs$gene_name[ix]

spe <- spe_list$mouseEmbryo_HVGs

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
g <- buildSNNGraph(spe, k = 50, use.dimred = "PCA")
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
ix <- which(res_list$mouseEmbryo_nnSVG$rank_nnSVG <= 1000)
top <- res_list$mouseEmbryo_nnSVG$gene_id[ix]

spe <- spe_list$mouseEmbryo_nnSVG

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
g <- buildSNNGraph(spe, k = 40, use.dimred = "PCA")
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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/res_downstream_clustering.rds")
saveRDS(res_out, file = fn)




