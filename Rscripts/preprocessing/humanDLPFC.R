######################
# Visium human DLPFC #
#     Load Data &    #
#    Preprocessing   #
######################


# Packages

library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scater)
library(scran)
library(here)


# ---------
# load data
# ---------

# load dataset as SpatialExperiment object from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)


# -------------
# preprocessing
# -------------

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)


# filter low-expressed and mitochondrial genes
# using default filtering parameters
spe <- filter_genes(spe)

dim(spe)


# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("/projectnb/weber-lr/lerusnak/outputs", "spe_humanDLPFC_preprocessed.rds")
saveRDS(spe, file = fn)