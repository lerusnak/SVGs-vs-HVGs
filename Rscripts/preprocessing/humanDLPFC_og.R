# Visium_Human_DLPFC

# Visium_Human_DLPFC

## Packages
library(SpatialExperiment)
library(STexampleData)
library(scater)
library(scran)
library(here)


## Load Data
humanDLPFC <- Visium_humanDLPFC()

# Pre-processing

## QC

###  save as spe
spe <- humanDLPFC

### library(scater)

### subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

### identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)

### calculate per-spot QC metrics and store in colData
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))

### select QC threshold for library size
qc_lib_size <- colData(spe)$sum < 600
colData(spe)$qc_lib_size <- qc_lib_size

### select QC threshold for number of expressed genes
qc_detected <- colData(spe)$detected < 400
colData(spe)$qc_detected <- qc_detected

### select QC threshold for mitochondrial read proportion
qc_mito <- colData(spe)$subsets_mito_percent > 28
colData(spe)$qc_mito <- qc_mito

### select QC threshold for number of cells per spot
qc_cell_count <- colData(spe)$cell_count > 10
colData(spe)$qc_cell_count <- qc_cell_count

### combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
### store in object
colData(spe)$discard <- discard

### remove combined set of low-quality spots
spe <- spe[, !colData(spe)$discard]


# filter zero-expressed genes
is_zero <- rowSums(counts(spe)) == 0
table(is_zero)

spe <- spe[!is_zero, ]

dim(spe)


## Normalization

### library(scran)

# calculate library size factors
spe <- computeLibraryFactors(spe)

# log-normalize
spe <-  logNormCounts(spe)

# remove mitochondrial genes
spe <- spe[!is_mito, ]


# -----------
# save object
# -----------

fn <- here("/projectnb/weber-lr/lerusnak/outputs", "spe_humanDLPFCog1.rds")
saveRDS(spe, file = fn)