###################
#    humanCRC     #
#   clustering    #
#      plots      #
###################


##############
#  packages  #
##############

library(SpatialExperiment)
library(here)
library(tidyverse)
library(mclust)
library(ggspavis)


###############
#  load data  #
###############

# directories
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_kmeansclust.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)



##################
# 'ground truth' #
#     labels     #
##################

## Create Results Dataframe

HVGs_df <- as.data.frame(cbind(coldata_out$HVGs, spatialcoords_out$HVGs)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
head(HVGs_df)

SPARKX_df <- as.data.frame(cbind(coldata_out$SPARKX, spatialcoords_out$SPARKX)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
head(SPARKX_df)

VGdf_list <- list(HVGs_df, SPARKX_df)



##################
# match clusters #
##################

# match clusters to ground truth layers for each method

#match_HVGs <- c(1, 3, 5, 4, 2)
#coldata_out$HVGs$label <- factor(
#  coldata_out$HVGs$label, levels = match_HVGs)

#match_SPARKX <- c(5, 1, 4, 3, 2)
#coldata_out$SPARKX$label <- factor(
#  coldata_out$SPARKX$label, levels = match_SPARKX)



###################
#  spatial plots  #
###################

VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- combined_df
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label)) + 
    geom_point(shape = ".", alpha = 0.5) + 
    coord_fixed() + 
    scale_y_reverse() + 
    #scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human CRC Layers: ", names(coldata_out)[i])) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanCRC_layers_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 6, height = 6.25)
}
