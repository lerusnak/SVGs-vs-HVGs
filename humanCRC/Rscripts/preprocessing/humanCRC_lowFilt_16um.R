##################
#    human CRC   #
# pre-processing #
##################


## Packages
library(SpatialExperiment)
library(scater)
library(scran)
library(here)
library(ggspavis)


#############
# Load Data #
#############

spe <- readRDS(file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_16um.rds")
dim(spe)

###############
# Filter Data #
###############

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)


#-------------------------------------------#
# NOTE: NO SPOT-BASED QC FILTERING OF SPOTS #
#-------------------------------------------#


# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rownames(spe))
table(is_mito)

# filter mitochondrial genes
spe <- spe[!is_mito, ]

dim(spe)


# filter low-expressed genes
is_low <- rowSums(counts(spe)) <= 10000
table(is_low)

spe <- spe[!is_low,]

dim(spe)

spe <- spe[, colSums(counts(spe)) > 0]

# QC 

spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
head(colData(spe))

hist(colData(spe)$sum, breaks = 20)


# Normalization 

# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)

assay(spe)[1:10, 1:10]
table(colData(spe)$sizeFactor >= 0)
anyMissing(colData(spe)$sizeFactor)

spe <- logNormCounts(spe)

assayNames(spe)

spe


#####################
# Save as .rds data #
#####################

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_lowFilt.rds")
saveRDS(spe, file = fn)




