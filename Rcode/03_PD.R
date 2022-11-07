
library(tidyverse)
library(rtrees)

#devtools::install_github("daijiang/lirrr")

#getting names
raptor_spp = colnames(comm)# species with distributions >1% of a pixel

# working with 200 trees downloaded from birdtree.org with Hackett backbone
# adding species that are not in phylogeny


raptors_trees = readRDS("Data/raw/raptors_trees.rds")

## 43 species added; 68 species added at genus level, 19 spp added at family level ##

#calculate Faiths PD
if (!file.exists("Data/output/pd_raptors.rds")) {
  n_phy = readRDS("Data/raw/raptors_trees.rds")
  comm_sub4 = as.data.frame(as.matrix(comm))#552 spp
  
  pd_list = purrr::map(n_phy, .f = function(x){
   lirrr::get_pd_alpha(samp_wide = comm_sub4, tree = x)
  }) # will take about 1 hour to run
  names(pd_list) = paste0("tree_", 1:200)
  pd_rapt = plyr::ldply(pd_list)
  pd_rapt = group_by(pd_rapt, site) %>% 
    summarise_if(is.numeric, mean, na.rm = T)
  saveRDS(pd_rapt, "Data/output/pd_raptors.rds")
} else {
  pd_rapt = readRDS("Data/output/pd_raptors.rds")
}

pd_rapt = pd_rapt %>%
  tidyr::separate(site, c("x", "y"), sep = "_") 




