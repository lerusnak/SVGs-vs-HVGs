######################
#   Human CRC Data   #
######################

# Packages
library(here)
library(arrow)
library(SpatialExperiment)
library(tibble)

#############################################

# data directory
CRC.outs <- "/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/data/outs"

# read in .parquet spatial coordinates file and convert to .csv
tissue_positions <- read_parquet(here(CRC.outs, "spatial", "tissue_positions.parquet"))
write.csv(tissue_positions, file = here(CRC.outs, "spatial", "tissue_positions_list.csv"))


# read in counts
fnm <- file.path(CRC.outs, "filtered_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)

# read in image data
img <- readImgData(path = file.path(CRC.outs, "spatial"), sample_id="crc")

# read in spatial coordinates
fnm <- file.path(CRC.outs, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm)

# construct observation & feature metadata
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)

# check 
length(intersect(colData(sce)$Barcode, xyz$barcode))
sum(duplicated(xyz$barcode))

# filter for spatial coords only in column data
dim(xyz)
xyz <- xyz[xyz$barcode %in% colData(sce)$Barcode,]
dim(xyz)

# add rownames to xyz
rownames(xyz) <- xyz$barcode
head(xyz)
tail(xyz)

# construct 'SpatialExperiment'
(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  rowData = rd,
  colData = DataFrame(xyz),
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img,
  sample_id = "crc"))

spe

# save spe as .rds file
saveRDS(spe, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC.rds")


