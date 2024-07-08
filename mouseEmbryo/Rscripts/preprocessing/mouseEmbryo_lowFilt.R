##################
#    human GBM   #
# pre-processing #
##################


## Packages
library(SpatialExperiment)
library(STexampleData)
library(scater)
library(scran)
library(here)


#############
# Load Data #
#############

spe <- readRDS(file = "/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/SlideSeqV2_mouseEmbryo.rds")
dim(spe)

###############
# Filter Data #
###############

# keep only spots over tissue - NO "IN_TISSUE" IN THIS DATASET
# spe <- spe[, colData(spe)$in_tissue == 1]
# dim(spe)


#-------------------------------------------#
# NOTE: NO SPOT-BASED QC FILTERING OF SPOTS #
#-------------------------------------------#


# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rownames(spe))
table(is_mito)

# filter mitochondrial genes
spe <- spe[!is_mito, ]

dim(spe)


# filter zero-expressed genes
# is_zero <- rowSums(counts(spe)) == 0
# spe <- spe[!is_zero, ]

# filter low-expressed genes
is_low <- rowSums(counts(spe)) <= 200 
table(is_low)

spe <- spe[!is_low,]

dim(spe)

# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

assayNames(spe)

spe


#####################
# Save as .rds data #
#####################

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_lowFilt.rds")
saveRDS(spe, file = fn)




