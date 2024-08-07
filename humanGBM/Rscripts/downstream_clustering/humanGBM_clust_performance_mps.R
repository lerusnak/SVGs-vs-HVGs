###################
#    humanGBM     #
#   clustering    #
#      plots      #
#  META PROGRAMS  #
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
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_clustering_mps.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)



##################
# 'ground truth' #
#     labels     #
##################

visium_metadata <- read.csv("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/visium_metadata.csv")

spot_labels <- visium_metadata %>%
  select(spot_id, sample, mp, layer) %>%
  filter(sample == "MGH258") 

dim(spot_labels)


## Create Results Dataframe

HVGs_df <- as.data.frame(cbind(coldata_out$HVGs, spatialcoords_out$HVGs)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
HVGs_df <- full_join(HVGs_df, spot_labels, by = "spot_id")

nnSVG_df <- as.data.frame(cbind(coldata_out$nnSVG, spatialcoords_out$nnSVG)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
nnSVG_df <- full_join(nnSVG_df, spot_labels, by = "spot_id")

SPARKX_df <- as.data.frame(cbind(coldata_out$SPARKX, spatialcoords_out$SPARKX)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
SPARKX_df <- full_join(SPARKX_df, spot_labels, by = "spot_id")

MorI_df <- as.data.frame(cbind(coldata_out$MorI, spatialcoords_out$MorI)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
MorI_df <- full_join(MorI_df, spot_labels, by = "spot_id")

VGdf_list <- list(HVGs_df, nnSVG_df, SPARKX_df, MorI_df)



##################
# match clusters #
##################

# match clusters to ground truth layers for each method

match_HVGs <- c(1, 2, 3, 6, 5, 4)
coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, levels = match_HVGs)

match_SPARKX <- c(1, 5, 3, 2, 4)
coldata_out$SPARKX$label <- factor(
  coldata_out$SPARKX$label, levels = match_SPARKX)

match_MorI <- c(3, 1, 2, 4, 5)
coldata_out$MorI$label <- factor(
  coldata_out$MorI$label, levels = match_MorI)

match_SDE2 <- c(2, 1, 3, 4, 5)
coldata_out$SpatialDE2$label <- factor(
  coldata_out$SpatialDE2$label, levels = match_SDE2)


###################
#  spatial plots  #
###################

VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- full_join(combined_df, spot_labels, by = "spot_id")
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label)) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    #scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human GBM (MPs): ", names(coldata_out)[i])) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanGBM_mps_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 4.25)
}



## Ground Truth plot

# plot ground truth MP labels in spatial coordinates
xyplot_MPs <-   ggplot(VGs_df[[1]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = mp)) + 
  geom_point(size = 1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  #scale_color_manual(values = pal, name = "label") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  ggtitle("Clustering", subtitle = "Human GBM: Meta Programs") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

ggsave("xyplot_trueMPs.png", path = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering")



###########
#   ARI   #
###########

## MPs ##

# calculate adjusted rand index

ari_HVGs_mps <- round(adjustedRandIndex(VGs_df[[1]]$label, 
                                  VGs_df[[1]]$mp), 4)

ari_nnSVG_mps <- round(adjustedRandIndex(VGs_df[[2]]$label, 
                                   VGs_df[[2]]$mp), 4)

ari_SPARKX_mps <-round(adjustedRandIndex(VGs_df[[3]]$label, 
                                            VGs_df[[3]]$mp), 4)

ari_MorI_mps <-round(adjustedRandIndex(VGs_df[[4]]$label, 
                                          VGs_df[[4]]$mp), 4)

ari_SDE2_mps <-round(adjustedRandIndex(VGs_df[[5]]$label, 
                                         VGs_df[[5]]$mp), 4)


# plot adjusted Rand index

df_mps <- data.frame(
  method = c("HVGs", "nnSVG", "SPARKX", "MoransI", "SpatialDE2"), 
  ARI = c(ari_HVGs_mps, ari_nnSVG_mps, ari_SPARKX_mps, ari_MorI_mps, ari_SDE2_mps)
)

df_mps


write_tsv(df_mps, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering/ari_mps.tsv")

pal_methods <- c("#D1A546", "#C459A1", "#50AD95", "#7D68FA", "#1570AD")

ggplot(df_mps, aes(x = method, y = ARI, shape = method, color = method)) + 
  geom_point(stroke = 1.5, size = 2) + 
  scale_shape_manual(values = c(4, 3, 5, 6, 1)) + 
  scale_color_manual(values = pal_methods) + 
  ylim(c(0, 1)) + 
  ggtitle("Downstream clustering performance") + 
  labs(y = "Adjusted Rand Index") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())

fn <- file.path(plots_dir, "summary_clustering_performance_ARI_mps")
ggsave(paste0(fn, ".png"), width = 4.75, height = 3.1)


