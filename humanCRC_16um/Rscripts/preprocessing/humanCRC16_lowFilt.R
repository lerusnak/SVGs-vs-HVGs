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

spe <- readRDS(file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs/spe_humanCRC16.rds")
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
is_low <- rowSums(counts(spe)) <= 100
table(is_low)

spe <- spe[!is_low,]

dim(spe)

# removing bins with zero counts
zero_bins <- colSums(counts(spe)) <= 0
table(zero_bins)

spe <- spe[, !zero_bins]
dim(spe)

# Normalization 

# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)

spe <- logNormCounts(spe)

assayNames(spe)

spe


#####################
# Save as .rds data #
#####################

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs/spe_humanCRC16_lowFilt.rds")
saveRDS(spe, file = fn)

