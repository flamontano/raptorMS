## FD ##

## MORPHOLOGICAL TRAITS 

row.names(morph_traits) = morph_traits$sp_map
morph_traits$sp_map = NULL

# community data with morphological information
comm_sub2 = comm[, colnames(comm) %in% row.names(morph_traits)]

# convert to presence/absence data 
comm_sub2[comm_sub2 < 1] = 0
comm_sub2[comm_sub2 >= 1] = 1
comm_sub2 = comm_sub2[rowSums(comm_sub2) > 0, ] # remove empty cells
comm_sub2 = comm_sub2[, colSums(comm_sub2) > 0] # remove empty species

morph_raptors = morph_traits[row.names(morph_traits) %in% colnames(comm_sub2),] # remove spp do not have distributional data

morph_raptors = scale(morph_raptors, center = T, scale = T)

if(!file.exists("Data/output/raptor_morph_FM.rds")){
  raptor_morph_fd = FD::dbFD(morph_raptors, as.matrix(comm_sub2), calc.FRic = F, 
                            calc.CWM = F, calc.FGR = F, calc.FDiv = F)
  saveRDS(raptor_morph_fd, file = "Data/output/raptor_morph_FM.rds")
} else {
  raptor_morph_fd = readRDS("Data/output/raptor_morph_FM.rds")
}


# null models ----
if(!file.exists("Data/output/morph_FM_null_df.rds")){
  morph_raptor_dist = dist(morph_raptors)
  morph_raptor_fdisp <- FD::fdisp(morph_raptor_dist, as.matrix(comm_sub2))
  
  ### 
  
  rap_morph_fdisp_null = parallel::mclapply(1:1000, function(i){
    set.seed(i)
    morph_raptors2 = morph_raptors
    row.names(morph_raptors2) = sample(row.names(morph_raptors2))
    # need to be the same order as comm_sub2
    morph_raptors2 = morph_raptors2[colnames(comm_sub2), ]
    morph_rap_dist2 = dist(morph_raptors2)
    FD::fdisp(morph_rap_dist2, as.matrix(comm_sub2))
  }, mc.cores = 55)
  
  rap_morph_fdisp_null_df = plyr::ldply(rap_morph_fdisp_null, function(x) x[[1]])
  rap_morph_fdisp_null_df = as.data.frame(t(rap_morph_fdisp_null_df))
  rap_morph_fdisp_null_df$fdis_obs = morph_raptor_fdisp$FDis
  
  # ut coord back
  rap_morph_fdisp_null_df$xy = row.names(rap_morph_fdisp_null_df)
  rap_morph_fdisp_null_df = tidyr::separate(rap_morph_fdisp_null_df, xy, c("x", "y"), sep = "_")
  
  # null mean and sd
  rap_morph_fdisp_null_df$null_mean = select(rap_morph_fdisp_null_df, starts_with("V")) %>% rowMeans() 
  rap_morph_fdisp_null_df$null_sd = apply(select(rap_morph_fdisp_null_df, starts_with("V")), 1, sd)
  
  all(rap_morph_fdisp_null_df$x %in% as.character(coord_continents$x))
  mean(rap_morph_fdisp_null_df$x %in% coord_continents$x)
  
  # merge continent info
  rap_morph_fdisp_null_df_terr = dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                                                    mutate(x = as.character(x), y = as.character(y)),
                                                  select(rap_morph_fdisp_null_df, x, y, fdis_obs, null_mean, null_sd),
                                                  by = c("x", "y")) %>% 
    mutate(x = as.numeric(x), y = as.numeric(y),
           fdis_ses = ifelse(null_sd > 0, (fdis_obs - null_mean)/null_sd, NA))
  # save summarised file
  saveRDS(rap_morph_fdisp_null_df_terr, file = "Data/output/morph_FM_null_df.rds")
} else {
  rap_morph_fdisp_null_df_terr = readRDS("Data/output/morph_FM_null_df.rds")
}


## NICHE TRAITS


row.names(niche_traits) = niche_traits$sp_map
niche_traits$sp_map = NULL

# community data with niche information
comm_sub2 = comm[, colnames(comm) %in% row.names(niche_traits)]


# convert to presence/absence data -- this deletes from 557 to 536 species . Restricted distributions?
comm_sub2[comm_sub2 < 1] = 0
comm_sub2[comm_sub2 >= 1] = 1
comm_sub2 = comm_sub2[rowSums(comm_sub2) > 0, ] # remove empty cells
comm_sub2 = comm_sub2[, colSums(comm_sub2) > 0] # remove empty species

niche_raptors = niche_traits[row.names(niche_traits) %in% colnames(comm_sub2),] # 

if(!file.exists("Data/output/raptor_niche_FM.rds")){
  niche_raptor_dist = FD::gowdis(niche_raptors)
  raptor_niche_fd <-  FD::dbFD(niche_raptor_dist, as.matrix(comm_sub2), calc.FRic = F, 
                             calc.CWM = F, calc.FGR = F, calc.FDiv = F)
  saveRDS(raptor_niche_fd, file = "Data/output/raptor_niche_FM.rds")
} else {
  raptor_niche_fd = readRDS("Data/output/raptor_niche_FM.rds")
}


# null models ----
if(!file.exists("Data/output/niche_FM_null_df.rds")){
  niche_raptor_dist = FD::gowdis(niche_raptors)
  niche_raptor_fdisp <- FD::fdisp(niche_raptor_dist, as.matrix(comm_sub2))
  
   rap_niche_fdisp_null = parallel::mclapply(1:1000, function(i){
    set.seed(i)
    niche_raptors2 = niche_raptors
    row.names(niche_raptors2) = sample(row.names(niche_raptors2))
    # need to be the same order as comm_sub2
    niche_raptors2 = niche_raptors2[colnames(comm_sub2), ]
    niche_rap_dist2 = FD::gowdis(niche_raptors2)
    FD::fdisp(niche_rap_dist2, as.matrix(comm_sub2))
  }, mc.cores = 55)
  
  rap_niche_fdisp_null_df = plyr::ldply(rap_niche_fdisp_null, function(x) x[[1]])
  rap_niche_fdisp_null_df = as.data.frame(t(rap_niche_fdisp_null_df))
  rap_niche_fdisp_null_df$fdis_obs = niche_raptor_fdisp$FDis
  
  # ut coord back
  rap_niche_fdisp_null_df$xy = row.names(rap_niche_fdisp_null_df)
  rap_niche_fdisp_null_df = tidyr::separate(rap_niche_fdisp_null_df, xy, c("x", "y"), sep = "_")
  
  # null mean and sd
  rap_niche_fdisp_null_df$null_mean = select(rap_niche_fdisp_null_df, starts_with("V")) %>% rowMeans() 
  rap_niche_fdisp_null_df$null_sd = apply(select(rap_niche_fdisp_null_df, starts_with("V")), 1, sd)
  
  all(rap_niche_fdisp_null_df$x %in% as.character(coord_continents$x))
  mean(rap_niche_fdisp_null_df$x %in% coord_continents$x)
  
  # merge continent info
  rap_niche_fdisp_null_df_terr2 = dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                                                     mutate(x = as.character(x), y = as.character(y)),
                                                   select(rap_niche_fdisp_null_df, x, y, fdis_obs, null_mean, null_sd),
                                                   by = c("x", "y")) %>% 
    mutate(x = as.numeric(x), y = as.numeric(y),
           fdis_ses = ifelse(null_sd > 0, (fdis_obs - null_mean)/null_sd, NA))
  # save summarized file
  saveRDS(rap_niche_fdisp_null_df_terr2, file = "Data/output/niche_FM_null_df.rds")
} else {
  rap_niche_fdisp_null_df_terr2 = readRDS("Data/output/niche_FM_null_df.rds")
}


## DISPERSAL TRAITS

row.names(dispersal_traits) = dispersal_traits$sp_map
dispersal_traits$sp_map = NULL

# community data with vagility information
comm_sub2 = comm[, colnames(comm) %in% row.names(dispersal_traits)]

# convert to presence/absence data -- this deletes from 557 to 536 species . Restricted distributions?
comm_sub2[comm_sub2 < 1] = 0
comm_sub2[comm_sub2 >= 1] = 1
comm_sub2 = comm_sub2[rowSums(comm_sub2) > 0, ] # remove empty cells
comm_sub2 = comm_sub2[, colSums(comm_sub2) > 0] # remove empty species

disp_raptors = dispersal_traits[row.names(dispersal_traits) %in% colnames(comm_sub2),] # 547

if(!file.exists("Data/output/raptor_dispersal_FM.rds")){
  disp_raptor_dist = FD::gowdis(disp_raptors)
  raptor_dispersal_fd <-  FD::dbFD(disp_raptor_dist, as.matrix(comm_sub2), calc.FRic = F, 
                               calc.CWM = F, calc.FGR = F, calc.FDiv = F)
  saveRDS(raptor_dispersal_fd , file = "Data/output/raptor_dispersal_FM.rds")
} else {
  raptor_dispersal_fd  = readRDS("Data/output/raptor_dispersal_FM.rds")
}


# null models ----
if(!file.exists("Data/output/dispersal_FM_null_df.rds")){
  disp_raptor_dist = FD::gowdis(disp_raptors)
  disp_raptor_fdisp <- FD::fdisp(disp_raptor_dist, as.matrix(comm_sub2))
  
  rap_disp_fdisp_null = parallel::mclapply(1:1000, function(i){
    set.seed(i)
    disp_raptors2 = disp_raptors
    row.names(disp_raptors2) = sample(row.names(disp_raptors2))
    # need to be the same order as comm_sub2
    disp_raptors2 = disp_raptors2[colnames(comm_sub2), ]
    disp_rap_dist2 = FD::gowdis(disp_raptors2)
    FD::fdisp(disp_rap_dist2, as.matrix(comm_sub2))
  }, mc.cores = 55)
  
  rap_disp_fdisp_null_df = plyr::ldply(rap_disp_fdisp_null, function(x) x[[1]])
  rap_disp_fdisp_null_df = as.data.frame(t(rap_disph_fdisp_null_df))
  rap_disp_fdisp_null_df$fdis_obs = disp_raptor_fdisp$FDis
  
  # ut coord back
  rap_disp_fdisp_null_df$xy = row.names(rap_disp_fdisp_null_df)
  rap_disp_fdisp_null_df = tidyr::separate(rap_disp_fdisp_null_df, xy, c("x", "y"), sep = "_")
  
  # null mean and sd
  rap_disp_fdisp_null_df$null_mean = select(rap_disp_fdisp_null_df, starts_with("V")) %>% rowMeans() 
  rap_disp_fdisp_null_df$null_sd = apply(select(rap_disp_fdisp_null_df, starts_with("V")), 1, sd)
  
  all(rap_disp_fdisp_null_df$x %in% as.character(coord_continents$x))
  mean(rap_disp_fdisp_null_df$x %in% coord_continents$x)
  
  # merge continent info
  rap_disp_fdisp_null_df_terr3 = dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                                                     mutate(x = as.character(x), y = as.character(y)),
                                                   select(rap_disp_fdisp_null_df, x, y, fdis_obs, null_mean, null_sd),
                                                   by = c("x", "y")) %>% 
    mutate(x = as.numeric(x), y = as.numeric(y),
           fdis_ses = ifelse(null_sd > 0, (fdis_obs - null_mean)/null_sd, NA))
  # save summarized file
  saveRDS(rap_disp_fdisp_null_df_terr3, file = "Data/output/dispersal_FM_null_df.rds")
} else {
  rap_disp_fdisp_null_df_terr3 = readRDS("Data/output/dispersal_FM_null_df.rds")
}

## DIET TRAITS ##

row.names(diet_traits) = diet_traits$sp_map
diet_traits$sp_map = NULL

# community data with diet information
comm_sub2 = comm[, colnames(comm) %in% row.names(diet_traits)]

# convert to presence/absence data -- this deletes from 557 to 536 species . Restricted distributions?
comm_sub2[comm_sub2 < 1] = 0
comm_sub2[comm_sub2 >= 1] = 1
comm_sub2 = comm_sub2[rowSums(comm_sub2) > 0, ] # remove empty cells
comm_sub2 = comm_sub2[, colSums(comm_sub2) > 0] # remove empty species

diet_raptors = diet_traits[row.names(diet_traits) %in% colnames(comm_sub2),] # 547

if(!file.exists("Data/output/raptor_diet_FM.rds")){
  diet_raptor_dist = FD::gowdis(diet_raptors)
  raptor_diet_fd <-  FD::dbFD(diet_raptor_dist, as.matrix(comm_sub2), calc.FRic = F, 
                                   calc.CWM = F, calc.FGR = F, calc.FDiv = F)
  saveRDS(raptor_diet_fd , file = "Data/output/raptor_diet_FM.rds")
} else {
  raptor_diet_fd  = readRDS("Data/output/raptor_diet_FM.rds")
}


# null models ----
if(!file.exists("Data/output/diet_FM_null_df.rds")){
  diet_raptor_dist = FD::gowdis(diet_raptors)
  diet_raptor_fdisp <- FD::fdisp(diet_raptor_dist, as.matrix(comm_sub2))
  
  rap_diet_fdisp_null = parallel::mclapply(1:1000, function(i){
    set.seed(i)
    diet_raptors2 = diet_raptors
    row.names(diet_raptors2) = sample(row.names(diet_raptors2))
    # need to be the same order as comm_sub2
    diet_raptors2 = diet_raptors2[colnames(comm_sub2), ]
    diet_rap_dist2 = FD::gowdis(diet_raptors2)
    FD::fdisp(diet_rap_dist2, as.matrix(comm_sub2))
  }, mc.cores = 55)
  
  rap_diet_fdisp_null_df = plyr::ldply(rap_diet_fdisp_null, function(x) x[[1]])
  rap_diet_fdisp_null_df = as.data.frame(t(rap_diet_fdisp_null_df))
  rap_diet_fdisp_null_df$fdis_obs = diet_raptor_fdisp$FDis
  
  # ut coord back
  rap_diet_fdisp_null_df$xy = row.names(rap_diet_fdisp_null_df)
  rap_diet_fdisp_null_df = tidyr::separate(rap_diet_fdisp_null_df, xy, c("x", "y"), sep = "_")
  
  # null mean and sd
  rap_diet_fdisp_null_df$null_mean = select(rap_diet_fdisp_null_df, starts_with("V")) %>% rowMeans() 
  rap_diet_fdisp_null_df$null_sd = apply(select(rap_diet_fdisp_null_df, starts_with("V")), 1, sd)
  
  all(rap_diet_fdisp_null_df$x %in% as.character(coord_continents$x))
  mean(rap_diet_fdisp_null_df$x %in% coord_continents$x)
  
  # merge continent info
  rap_diet_fdisp_null_df_terr4 = dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                                                      mutate(x = as.character(x), y = as.character(y)),
                                                    select(rap_diet_fdisp_null_df, x, y, fdis_obs, null_mean, null_sd),
                                                    by = c("x", "y")) %>% 
    mutate(x = as.numeric(x), y = as.numeric(y),
           fdis_ses = ifelse(null_sd > 0, (fdis_obs - null_mean)/null_sd, NA))
  # save summarized file
  saveRDS(rap_diet_fdisp_null_df_terr4, file = "Data/output/diet_FM_null_df.rds")
} else {
  rap_diet_fdisp_null_df_terr4 = readRDS("Data/output/diet_FM_null_df.rds")
}

## FORAGING TRAITS ##

row.names(forag_traits) = forag_traits$sp_map
forag_traits$sp_map = NULL

# community data with foraging information
comm_sub2 = comm[, colnames(comm) %in% row.names(forag_traits)]


# convert to presence/absence data
comm_sub2[comm_sub2 < 1] = 0
comm_sub2[comm_sub2 >= 1] = 1
comm_sub2 = comm_sub2[rowSums(comm_sub2) > 0, ] # remove empty cells
comm_sub2 = comm_sub2[, colSums(comm_sub2) > 0] # remove empty species

forag_raptors = forag_traits[row.names(forag_traits) %in% colnames(comm_sub2),] # 547

if(!file.exists("Data/output/raptor_forag_FM.rds")){
  forag_raptor_dist = FD::gowdis(forag_raptors)
  raptor_forag_fd <-  FD::dbFD(forag_raptor_dist, as.matrix(comm_sub2), calc.FRic = F, 
                              calc.CWM = F, calc.FGR = F, calc.FDiv = F)
  saveRDS(raptor_forag_fd , file = "Data/output/raptor_forag_FM.rds")
} else {
  raptor_forag_fd  = readRDS("Data/output/raptor_forag_FM.rds")
}


# null models ----
if(!file.exists("Data/output/forag_FM_null_df.rds")){
  forag_raptor_dist = FD::gowdis(forag_raptors)
  forag_raptor_fdisp <- FD::fdisp(forag_raptor_dist, as.matrix(comm_sub2))
  
  rap_forag_fdisp_null = parallel::mclapply(1:1000, function(i){
    set.seed(i)
    forag_raptors2 = forag_raptors
    row.names(forag_raptors2) = sample(row.names(forag_raptors2))
    # need to be the same order as comm_sub2
    forag_raptors2 = forag_raptors2[colnames(comm_sub2), ]
    forag_rap_dist2 = FD::gowdis(forag_raptors2)
    FD::fdisp(morph_rap_dist2, as.matrix(comm_sub2))
  }, mc.cores = 55)
  
  rap_forag_fdisp_null_df = plyr::ldply(rap_forag_fdisp_null, function(x) x[[1]])
  rap_forag_fdisp_null_df = as.data.frame(t(rap_forag_fdisp_null_df))
  rap_forag_fdisp_null_df$fdis_obs = forag_raptor_fdisp$FDis
  
  # ut coord back
  rap_forag_fdisp_null_df$xy = row.names(rap_forag_fdisp_null_df)
  rap_forag_fdisp_null_df = tidyr::separate(rap_forag_fdisp_null_df, xy, c("x", "y"), sep = "_")
  
  # null mean and sd
  rap_forag_fdisp_null_df$null_mean = select(rap_forag_fdisp_null_df, starts_with("V")) %>% rowMeans() 
  rap_forag_fdisp_null_df$null_sd = apply(select(rap_forag_fdisp_null_df, starts_with("V")), 1, sd)
  
  all(rap_forag_fdisp_null_df$x %in% as.character(coord_continents$x))
  mean(rap_forag_fdisp_null_df$x %in% coord_continents$x)
  
  # merge continent info
  rap_forag_fdisp_null_df_terr5 = dplyr::right_join(dplyr::select(coord_continents, x, y, country, continent) %>% 
                                                      mutate(x = as.character(x), y = as.character(y)),
                                                    select(rap_forag_fdisp_null_df, x, y, fdis_obs, null_mean, null_sd),
                                                    by = c("x", "y")) %>% 
    mutate(x = as.numeric(x), y = as.numeric(y),
           fdis_ses = ifelse(null_sd > 0, (fdis_obs - null_mean)/null_sd, NA))
  # save summarized file
  saveRDS(rap_forag_fdisp_null_df_terr5, file = "Data/output/forag_FM_null_df.rds")
} else {
  rap_forag_fdisp_null_df_terr5 = readRDS("Data/output/forag_FM_null_df.rds")
}


# binding #

names(rap_morph_fdisp_null_df_terr)[5:8] = paste0(names(rap_morph_fdisp_null_df_terr)[5:8], "_morph")
names(rap_niche_fdisp_null_df_terr2)[5:8] = paste0(names(rap_niche_fdisp_null_df_terr2)[5:8], "_niche")
names(rap_disp_fdisp_null_df_terr3)[5:8] = paste0(names(rap_disp_fdisp_null_df_terr3)[5:8], "_dispe")
names(rap_diet_fdisp_null_df_terr4)[5:8] = paste0(names(rap_diet_fdisp_null_df_terr4)[5:8], "_diet")
names(rap_forag_fdisp_null_df_terr5)[5:8] = paste0(names(rap_forag_fdisp_null_df_terr5)[5:8], "_forag")

rap_fd0 = left_join(rap_morph_fdisp_null_df_terr, 
                   rap_niche_fdisp_null_df_terr2,
                        by = c("x", "y", "country", "continent"))

rap_fd1 = left_join(rap_fd0, 
                   rap_disp_fdisp_null_df_terr3,
                   by = c("x", "y", "country", "continent"))
rap_fd2 = left_join(rap_fd1, 
                   rap_diet_fdisp_null_df_terr4,
                   by = c("x", "y", "country", "continent"))
rap_fd3 = left_join(rap_fd2, 
                    rap_forag_fdisp_null_df_terr5,
                    by = c("x", "y", "country", "continent"))
rap_fd<-rap_fd3


names(rap_fd)

