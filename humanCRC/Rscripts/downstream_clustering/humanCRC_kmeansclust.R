############################
#        human CRC         #
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


# methods: SPARK-X, SpatialDE2, HVGs

spe_list <- list(
  humanCRC_HVGs = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_HVGs.rds")),
  humanCRC_SPARKX = readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_SPARKX.rds"))
  )

res_list <- list(
  humanCRC_SPARKX = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_SPARKX.rds"))),
  humanCRC_HVGs = rowData(readRDS(here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanCRC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanCRC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanCRC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanCRC_HVGs"]]), "_HVGs")[-(1:2)]


table(res_list$humanCRC_HVGs$gene_name %in% res_list$humanCRC_SPARKX$gene_name)


spe_out <- list()


# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanCRC_HVGs$rank_HVGs <= 1000)
top <- res_list$humanCRC_HVGs$gene_name[ix]

spe <- spe_list$humanCRC_HVGs


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_hvgs <- kmeans(pcs, centers = 3)

table(km_hvgs$cluster)

# kmeans cluster assignments
km_hvg_11clust <- km_hvgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_hvg_11clust)


# store object
spe_out$spe_HVGs <- spe



# ------------------------------
# downstream clustering: SPARK-X
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# SPARK-X: top 1000 SVGs
ix <- which(res_list$humanCRC_SPARKX$rank_SPARKX <= 1000)
top <- res_list$humanCRC_SPARKX$gene_name[ix]

spe <- spe_list$humanCRC_SPARKX


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)

# Vector of principal commponents for clustering 
pcs <- reducedDim(spe, "PCA")


##### K means Clustering


set.seed(123)

km_svgs <- kmeans(pcs, centers = 3)

table(km_svgs$cluster)

# kmeans cluster assignments
km_svg_11clust <- km_svgs$cluster

# store cluster label in column 'label' in colData
colLabels(spe) <- factor(km_svg_11clust)


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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/res_downstream_kmeansclust.rds")
saveRDS(res_out, file = fn)

