###################
#  mouse Embryo   #
#   clustering    #
#     plots       #
###################


##############
#  packages  #
##############

library(SpatialExperiment)
library(here)
library(mclust)
library(ggplot2)
library(ggspavis)
library(readr)


###############
#  load data  #
###############

# directories
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_clustering.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


##################
# match clusters #
##################

# match clusters to ground truth layers for each method

#match_HVGs <- c(2, 4, 6, 3, 8, 7, 1, 5)
#coldata_out$HVGs$label <- factor(
#  coldata_out$HVGs$label, levels = match_HVGs)

#match_nnSVG <- c(6, 7, 5, 8, 2, 4, 3, 1)
#coldata_out$nnSVG$label <- factor(
#  coldata_out$nnSVG$label, levels = match_nnSVG)


###################
#  spatial plots  #
###################

for (i in seq_along(coldata_out)) {
  
  df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]]))
  
  ggplot(df, aes(x = X2, y = X3, 
                 color = label)) + 
    geom_point(size = 0.1, alpha = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Mouse Embryo: ", names(coldata_out[i]))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(plots_dir, 
                  paste0("clustering_mouseEmbryo_", names(coldata_out[i])))
  ggsave(paste0(fn, ".pdf"), width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}



### Compositional Analysis

library(reshape)


length(coldata_out[["HVGs"]]$label)
length(coldata_out[["nnSVG"]]$label)

# check cluster frequencies
table(coldata_out[["HVGs"]]$label)
table(coldata_out[["nnSVG"]]$label)


# get cluster labels
hvgs_labels <- sort(unique(coldata_out[["HVGs"]]$label))
nnsvg_labels <- sort(unique(coldata_out[["nnSVG"]]$label))

hvgs_labels
nnsvg_labels



# create empty matrix to store results
res <- matrix(
  NA, 
  ncol = length(hvgs_labels), 
  nrow = length(nnsvg_labels)
)

colnames(res) <- paste0("hvgs_", hvgs_labels)
rownames(res) <- paste0("nnsvg_", nnsvg_labels)

# convert vectors of cluster labels to factors so empty levels are not dropped when subsetting
hvgs_clus <- as.factor(coldata_out[["HVGs"]]$label)
nnsvg_clus <- as.factor(coldata_out[["nnSVG"]]$label)

# calculate by row (i.e. for each nnsvg domain)
for (i in 1:nrow(res)) {
  # subset hvgs cluster labels for points in nnsvg domain i
  hvgs_sub <- hvgs_clus[nnsvg_clus == nnsvg_labels[i]]
  # calculate frequency table
  hvgs_sub_tbl <- table(hvgs_sub)
  # convert frequency table to proportions and store in matrix
  res[i, ] <- hvgs_sub_tbl / sum(hvgs_sub_tbl)
}

res
rowSums(res)


df <- melt(cbind(nnsvg_labels, as.data.frame(res)), id.vars = "nnsvg_labels")
df



# then create stacked bar plot
ggplot(df, aes(x=nnsvg_labels, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  labs(y="Proportion of HVGs Cluster", 
       x="nnsvg Cluster",
       title = "Proportion of hvgs cluster spots \nwithin a nnsvg region") +
  theme_bw()


################################################

## Proportions of beads assigned to each cluster comapering VG methods

# HVG Proportions
hvg_props <- as.data.frame(table(coldata_out[["HVGs"]]$label)/length(coldata_out[["HVGs"]]$label))
colnames(hvg_props) <- c("cluster", "proportion")
hvg_props$proportion <- round(hvg_props$proportion, 4)
hvg_props

method <- rep("HVGs", 11)

hvg_props <- cbind(hvg_props, method)
hvg_props

# SVG Proportions
nnsvg_props <- as.data.frame(table(coldata_out[["nnSVG"]]$label)/length(coldata_out[["nnSVG"]]$label))
colnames(nnsvg_props) <- c("cluster", "proportion")
nnsvg_props$proportion <- round(nnsvg_props$proportion, 4)
nnsvg_props

method <- rep("nnSVG", 11)

nnsvg_props <- cbind(nnsvg_props, method)
nnsvg_props

# Combined VG proportions data frame
vg_props <- rbind(hvg_props, nnsvg_props)
vg_props

vg_props_wide <- vg_props %>% pivot_wider(names_from = cluster, values_from = proportion) %>% 
  select(method, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`, `11`)
vg_props_wide

write_tsv(vg_props_wide, file = here(plots_dir, 'prop_spots_percluster.tsv'))


## Frequencies

# HVG Frequencies
hvg_freqs <- as.data.frame(table(coldata_out[["HVGs"]]$label))
colnames(hvg_freqs) <- c("cluster", "frequency")
hvg_freqs

method <- rep("HVGs", 11)

hvg_freqs <- cbind(hvg_freqs, method)
hvg_freqs

# SVG Frequencies
nnsvg_freqs <- as.data.frame(table(coldata_out[["nnSVG"]]$label))
colnames(nnsvg_freqs) <- c("cluster", "frequency")
nnsvg_freqs

method <- rep("nnSVG", 11)

nnsvg_freqs <- cbind(nnsvg_freqs, method)
nnsvg_freqs

# Combined VG frequencies data frame
vg_freqs <- rbind(hvg_freqs, nnsvg_freqs)
vg_freqs

vg_freqs_wide <- vg_freqs %>% pivot_wider(names_from = cluster, values_from = frequency) %>% 
  select(method, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`, `11`)
vg_freqs_wide

write_tsv(vg_freqs_wide, file = here(plots_dir, 'beadsfreq_percluster.tsv'))

# Stacked Barplot
ggplot(data = vg_freqs, aes(x = method, y = frequency, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "frequency of beads belonging to ")
ggsave(filename = here(plots_dir, 'beadsfreq_percluster_plot.png'))


##############################################################

# Comparing HVG cluster labels to nnSVG cluster labels

nnsvg_vs_hvg_ari <- adjustedRandIndex(coldata_out[["nnSVG"]]$label,coldata_out[["HVGs"]]$label)
nnsvg_vs_hvg_ari # 0.558479

