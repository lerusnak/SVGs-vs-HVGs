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
fn <- here(clustering_dir, 'res_downstream_mbkmeansclust.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


###################
#  spatial plots  #
###################

## 5 clusters

# match clusters to ground truth layers for each method

match_HVGs <- c(5, 4, 3, 2, 1)
coldata_out$HVGs$label_k5 <- factor(
  coldata_out$HVGs$label_k5, levels = match_HVGs)
match_SPARKX <- c(5, 4, 1, 3, 2)
coldata_out$SPARKX$label_k5 <- factor(
  coldata_out$SPARKX$label_k5, levels = match_SPARKX)


VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- combined_df
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label_k5)) + 
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
          axis.ticks = element_blank()) +
    labs(color = "cluster") +
    guides(colour = guide_legend(override.aes = list(shape = 19)))
  
  print(p)
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanCRC_05_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 6, height = 6.25)
}

ari_SPARKX_5 <- adjustedRandIndex(VGs_df[[1]]$label_k5, 
                                       VGs_df[[2]]$label_k5)
ari_SPARKX_5


##########################################

## 10 clusters

# match clusters to ground truth layers for each method

match_HVGs <- c(3,8,7,4,1,6,5,2,9,10)
coldata_out$HVGs$label_k10 <- factor(
  coldata_out$HVGs$label_k10, levels = match_HVGs)
match_SPARKX <- c(2,5,7,10,8,9,6,1,4,3)
coldata_out$SPARKX$label_k10 <- factor(
  coldata_out$SPARKX$label_k10, levels = match_SPARKX)


VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- combined_df
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label_k10)) + 
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
          axis.ticks = element_blank()) +
    labs(color = "cluster") +
    guides(colour = guide_legend(override.aes = list(shape = 19)))
  
  print(p)
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanCRC_10_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 6, height = 6.25)
}

ari_SPARKX_10 <- adjustedRandIndex(VGs_df[[1]]$label_k10, 
                                        VGs_df[[2]]$label_k10)
ari_SPARKX_10


##########################################

## 15 clusters

# match clusters to ground truth layers for each method

match_HVGs <- c(7,12,2,8,4,9,1,15,10,6,11,3,13,14,5)
coldata_out$HVGs$label_k15 <- factor(
  coldata_out$HVGs$label_k15, levels = match_HVGs)
match_SPARKX <- c(15,10,5,1,6,4,11,9,13,12,8,7,14,2,3)
coldata_out$SPARKX$label_k15 <- factor(
  coldata_out$SPARKX$label_k15, levels = match_SPARKX)

VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- combined_df
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label_k15)) + 
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
          axis.ticks = element_blank()) +
    labs(color = "cluster") +
    guides(colour = guide_legend(override.aes = list(shape = 19), ncol=2))
  
  print(p)
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanCRC_15_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 6, height = 6.25)
}

ari_SPARKX_15 <- adjustedRandIndex(VGs_df[[1]]$label_k15, 
                                        VGs_df[[2]]$label_k15)
ari_SPARKX_15


## 20 clusters

# match clusters to ground truth layers for each method

match_HVGs <- c(5,7,13,4,1,6,2,8,9,18,11,12,15,14,16,3,17,10,19,20)
coldata_out$HVGs$label_k20 <- factor(
  coldata_out$HVGs$label_k20, levels = match_HVGs)
match_SPARKX <- c(8,12,3,4,13,6,15,18,9,10,2,11,19,14,7,16,17,5,1,20)
coldata_out$SPARKX$label_k20 <- factor(
  coldata_out$SPARKX$label_k20, levels = match_SPARKX)

VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- combined_df
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label_k20)) + 
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
          axis.ticks = element_blank()) +
    labs(color = "cluster") +
    guides(colour = guide_legend(override.aes = list(shape = 19), ncol=2))
  
  print(p)
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanCRC_20_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 6, height = 6.25)
}

ari_SPARKX_20 <- adjustedRandIndex(VGs_df[[1]]$label_k20, 
                                        VGs_df[[2]]$label_k20)
ari_SPARKX_20

