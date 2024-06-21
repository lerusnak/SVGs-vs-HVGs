###############
#  human GBM  #
#  Load Data  #
###############


# Packages
library(Biobase)
library(Seurat)
library(SpatialExperiment)
library(S4Vectors)


# Load Data 

# load as a seurat object

gbm_mgh258 <- Load10X_Spatial(data.dir = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/outs",
                              filename = "filtered_feature_bc_matrix.h5",
                              assay = "Spatial",
                              slice = "detected_tissue_image.jpg",
                              filter.matrix = TRUE,
                              to.upper = FALSE)


# Function to
# coerce data into spatial experiment object from seurat object

seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}


#####################################

spe <- seurat_to_spe(gbm_mgh258, sample_id = "MGH258", img_id = 1)

spe

saveRDS(spe, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM.rds")