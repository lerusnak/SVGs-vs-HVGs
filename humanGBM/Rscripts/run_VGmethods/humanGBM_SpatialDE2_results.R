#############################
#         Human GBM         #
# Format SpatialDE2 Results #
#############################

# packages
library(here)
library(dplyr)


# Load results from SpatialDE2 method run
fn.spatialde2 <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/humanGBM_SpatialDE2_results.csv")
rowdat_spatialde2 <- read.csv(fn.spatialde2)

sde2_ranks <- rowdat_spatialde2 %>% 
  arrange(desc(FSV)) %>% 
  mutate(rank = c(1:nrow(rowdat_spatialde2)), gene_id = X)

head(sde2_ranks)

write.csv(sde2_ranks, "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/humanGBM_SpatialDE2.csv")
