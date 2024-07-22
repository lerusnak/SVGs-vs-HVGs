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
  mouseEmbryo_nnSVG = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_nnSVG.rds")),
  mouseEmbryo_SPARKX = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_SPARKX.rds")),
  mouseEmbryo_MorI = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_MoransI.rds"))
)

res_list <- list(
  mouseEmbryo_nnSVG = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_nnSVG.rds"))),
  mouseEmbryo_SPARKX = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_SPARKX.rds"))),
  mouseEmbryo_MorI = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_MoransI.rds"))),
  mouseEmbryo_HVGs = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_HVGs.rds")))
  )

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseEmbryo_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["mouseEmbryo_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["mouseEmbryo_MorI"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_MorI"]]), "_MorI")[-(1:2)]
colnames(res_list[["mouseEmbryo_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["mouseEmbryo_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method

table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)
table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_SPARKX$gene_name)
table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_MorI$gene_name)


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

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_hvgs <- kmeans(pcs, centers = 11)

table(km_hvgs$cluster)

# kmeans cluster assignments
km_hvg_11clust <- km_hvgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_hvg_11clust)


# store object
spe_out$spe_HVGs <- spe



# ----------------------------
# downstream clustering: nnSVG
# ----------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# nnSVG: top 1000 SVGs
ix <- which(res_list$mouseEmbryo_nnSVG$rank_nnSVG <= 1000)
top <- res_list$mouseEmbryo_nnSVG$gene_name[ix]

spe <- spe_list$mouseEmbryo_nnSVG


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_svgs <- kmeans(pcs, centers = 11)

table(km_svgs$cluster)

# kmeans cluster assignments
km_svg_11clust <- km_svgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_svg_11clust)


# store object
spe_out$spe_nnSVG <- spe




# ------------------------------
# downstream clustering: SPARK-X
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# SPARK-X: top 1000 SVGs
ix <- which(res_list$mouseEmbryo_SPARKX$rank_SPARKX <= 1000)
top <- res_list$mouseEmbryo_SPARKX$gene_name[ix]

spe <- spe_list$mouseEmbryo_SPARKX


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_svgs <- kmeans(pcs, centers = 11)

table(km_svgs$cluster)

# kmeans cluster assignments
km_svg_11clust <- km_svgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_svg_11clust)


# store object
spe_out$spe_SPARKX <- spe



# ----------------------------
# downstream clustering: MorI
# ----------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# Moran's I: top 1000 SVGs
ix <- which(res_list$mouseEmbryo_MorI$rank_MorI <= 1000)
top <- res_list$mouseEmbryo_MorI$gene_name[ix]

spe <- spe_list$mouseEmbryo_MorI


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_svgs <- kmeans(pcs, centers = 11)

table(km_svgs$cluster)

# kmeans cluster assignments
km_svg_11clust <- km_svgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_svg_11clust)


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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/res_downstream_kmeansclust.rds")
saveRDS(res_out, file = fn)




