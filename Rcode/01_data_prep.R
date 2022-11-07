
### DATA #### -------------------------------------------------------------------------------------------------------

#Importing data #

## trait data
traits<-read.csv("Data/raw/raptor_traits.csv", header=T)#557 species

#morphological
morph_traits<-traits[,c(1:9)];morph_traits<-morph_traits[complete.cases(morph_traits), ]#556 species

#ecological niche
niche_traits<-traits[,c(1,10,12)];niche_traits<-niche_traits[complete.cases(niche_traits), ]#556
  niche_traits$Primary.Lifestyle <-as.factor(niche_traits$Primary.Lifestyle)

  #diet
diet_traits<-traits[,c(1, 14:20)];diet_traits<-diet_traits[complete.cases(diet_traits), ]#468 

#summary(diet_traits)

# foraging
forag_traits<-traits[,c(1, 21:27)];forag_traits<-forag_traits[complete.cases(forag_traits), ]#468

# vagility
dispersal_traits<-traits[,c(1,11,13)];dispersal_traits<-dispersal_traits[complete.cases(dispersal_traits), ]#555
dispersal_traits$Migration<-as.numeric(dispersal_traits$Migration)

## community data

comm = readRDS("Data/raw/raptor_distribution_eck4_all_Matrix.rds")

# coordinates of cells

coord_continents = readRDS("Data/output/coord_continent.rds") %>% as_tibble()




