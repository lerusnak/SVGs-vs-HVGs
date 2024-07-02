######################
#     human GBM      #
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


humanGBM_HVGs <- humanGBM_HVGs[,colnames(humanGBM_HVGs) %in% spots]
dim(humanGBM_HVGs)

humanGBM_nnSVG <- humanGBM_nnSVG[,colnames(humanGBM_nnSVG) %in% spots]
dim(humanGBM_nnSVG)


spe_list <- list(humanGBM_HVGs = humanGBM_HVGs, humanGBM_nnSVG = humanGBM_nnSVG)

res_list <- list(humanGBM_HVGs = rowData(humanGBM_HVGs), humanGBM_nnSVG = rowData(humanGBM_nnSVG))


# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanGBM_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanGBM_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanGBM_HVGs"]]), "_HVGs")[-(1:2)]


table(res_list$humanGBM_HVGs$symbol %in% res_list$humanGBM_nnSVG$symbol)

spe_out <- list()



# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanGBM_HVGs$rank_HVGs <= 1000)
genes_subset <- res_list$humanGBM_HVGs[ix,]
top <- rownames(genes_subset)

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
genes_subset <- res_list$humanGBM_nnSVG[ix,]
top <- rownames(genes_subset)

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
g <- buildSNNGraph(spe, k = 8, use.dimred = "PCA")
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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/res_downstream_clustering_layers.rds")
saveRDS(res_out, file = fn)




