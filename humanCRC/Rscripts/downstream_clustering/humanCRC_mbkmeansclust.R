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
library(mbkmeans)

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


table(rownames(res_list$humanCRC_HVGs) %in% rownames(res_list$humanCRC_SPARKX))


spe_out <- list()


##################
#   clustering   #
##################

# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanCRC_HVGs$rank_HVGs <= 1000)
top <- rownames(res_list$humanCRC_HVGs)[ix]

spe <- spe_list$humanCRC_HVGs


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)


# kmeans clustering

library(mbkmeans)


## 5 clusters

res <- mbkmeans(spe, clusters = 5,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k5 <- factor(res$Clusters)
head(colData(spe))


## 10 clusters

res <- mbkmeans(spe, clusters = 10,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k10 <- factor(res$Clusters)
head(colData(spe))


## 15 clusters

res <- mbkmeans(spe, clusters = 15,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k15 <- factor(res$Clusters)
head(colData(spe))


## 20 clusters

res <- mbkmeans(spe, clusters = 20,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k20 <- factor(res$Clusters)
head(colData(spe))


# store object
spe_out$spe_HVGs <- spe


# -------------------------------
# downstream clustering: SPARK-X
# -------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 SVGs
ix <- which(res_list$humanCRC_SPARKX$rank_SPARKX <= 1000)
top <- rownames(res_list$humanCRC_SPARKX)[ix]

spe <- spe_list$humanCRC_SPARKX


# dimensionality reduction

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)


# kmeans clustering

library(mbkmeans)

## 5 clusters

res <- mbkmeans(spe, clusters = 5,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k5 <- factor(res$Clusters)
head(colData(spe))


## 10 clusters

res <- mbkmeans(spe, clusters = 10,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k10 <- factor(res$Clusters)
head(colData(spe))


## 15 clusters

res <- mbkmeans(spe, clusters = 15,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k15 <- factor(res$Clusters)
head(colData(spe))


## 20 clusters

res <- mbkmeans(spe, clusters = 20,
                reduceMethod = "PCA",
                whichAssay = "logcounts")
head(res$Clusters)

# store cluster label in column 'label' in colData
colData(spe)$label_k20 <- factor(res$Clusters)
head(colData(spe))

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


fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/res_downstream_mbkmeansclust.rds")
saveRDS(res_out, file = fn)

