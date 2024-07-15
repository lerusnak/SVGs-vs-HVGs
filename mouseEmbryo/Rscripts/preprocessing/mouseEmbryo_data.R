####################
#   Mouse Embryo   #
####################

library(SpatialExperiment)
library(readr)
library(dplyr)


# ---------
# Load data
# ---------

# load data object Sampath Kumar, A., Tian, L., Bolondi, A. et al. Spatiotemporal transcriptomic maps of whole mouse embryos at the onset of organogenesis. Nat Genet 55, 1176–1185 (2023). https://doi.org/10.1038/s41588-023-01435-6
# previously downloaded and saved locally from GEO

dir_data <- "/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/data"


# load expression matrix
# note runtime: several minutes
exprs <- read_table(file.path(dir_data, "GSM5915040_201104_07.digital_expression.txt.gz"))

dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]

# get gene IDs from first column
gene_ids <- exprs$GENE
str(gene_ids)
length(gene_ids)

# get bead IDs from column names (excluding gene IDs column)
bead_ids <- colnames(exprs)[-1]
str(bead_ids)
length(bead_ids)

# convert expression matrix to numeric matrix without gene IDs
exprs <- exprs[, -1]
exprs <- as.matrix(exprs)
stopifnot(nrow(exprs) == length(gene_ids))
rownames(exprs) <- gene_ids

dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]


# load beads info
# note: contains mix of tab-delimited and comma-delimited
# "GSM5915058_201104_33_matched_bead_locations.txt.gz" vs "Puck_200115_08_bead_locations.csv"

file <- file.path(dir_data, "GSM5915040_201104_07_matched_bead_locations.txt.gz")

bead_locations <- read_tsv(file, col_names = FALSE)


dim(bead_locations)
head(bead_locations)

stopifnot(nrow(bead_locations) == ncol(exprs))


# ------------------------
# Create SpatialExperiment
# ------------------------

# convert assay matrix to sparse format
exprs_sparse <- as(exprs, "dgCMatrix")
format(object.size(exprs_sparse), units = "GB")

# row data
row_data <- DataFrame(gene_name = gene_ids)
rownames(row_data) <- row_data$gene_name

# column data
col_data <- DataFrame(
  barcode_id = bead_ids, 
  sample_id = "sample01"
)
rownames(col_data) <- col_data$barcode_id


# spatial coordinates
spatial_coords <- as.matrix(bead_locations[, c("X2", "X3")])
rownames(spatial_coords) <- bead_ids


spe <- SpatialExperiment(
  assays = list(counts = exprs_sparse), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords
)

spe


# ----------------
# Save data object
# ----------------

saveRDS(spe, file = "/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/SlideSeqV2_mouseEmbryo.rds")


