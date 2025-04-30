# Load packages and functions
library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
source("CredentialsAndSecurity/Cryption.R")
# Load crypting key and credentials
load("CredentialsAndSecurity/key.RData")
credentials <- read.aes(filename = "CredentialsAndSecurity/credentials.txt",key = key)

# Download file from the cloud
# Should be run once or when the file is updated
# webdav::webdav_download_file(
#   "https://nextcloud.inrae.fr/remote.php/dav/files/tlafontrapn",  #next cloud webdav adress
#   file_path = "INHERSEED/monitoring_and_info.xlsx", #path to the file I want to dl
#   destination_path = "test2beignored/data", #path to the local folder
#   username = credentials$login, #credentials created beforehand
#   password = credentials$password,
#   verbose = TRUE
# )
# webdav::webdav_download_file(
#   "https://nextcloud.inrae.fr/remote.php/dav/files/tlafontrapn",  #next cloud webdav adress
#   file_path = "INHERSEED/phenotype.xlsx", #path to the file I want to dl
#   destination_path = "test2beignored/data", #path to the local folder
#   username = credentials$login, #credentials created beforehand
#   password = credentials$password,
#   verbose = TRUE
# )

# Read file
data_ordered_seeds <- readxl::read_xlsx("test2beignored/data/monitoring_and_info.xlsx",sheet = 2) #list of ordered seeds
data_quantity_seeds <- readxl::read_xlsx("test2beignored/data/monitoring_and_info.xlsx",sheet = 5) #seeds sowing and stock tracking
data_phenotype_plants <- readxl::read_xlsx("test2beignored/data/phenotype.xlsx",sheet = 1) # phenotype info

# Format taxonomy

data_ordered_seeds%<>% #create new taxonomy columns
  mutate(genus=stringr::str_extract(species_for_tree,"^\\S*"),
         species=stringr::str_extract(species_for_tree,"\\S*$"))
data_ordered_seeds %<>%
  mutate(species=ifelse(nchar(species)<1,"sp",species))%>%
  mutate(species=paste0(genus,"_",species))
data_ordered_seeds%<>% #convert to numeric, will generate NA because some are text but its temporary
  mutate(nb_seeds_lot=as.numeric(nb_seeds_lot),
         weight_1000_seeds_g=as.numeric(weight_1000_seeds_g))%>%
  arrange(desc(weight_1000_seeds_g))%>%
  mutate(species=forcats::fct_inorder(species))

# Seed mass
data_ordered_seeds%>%
  ggplot(aes(x=species,y=weight_1000_seeds_g))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  ylab("1000 seeds weight (g)")+
  data_ordered_seeds%>%
  mutate(nb_seeds_needed=0.5/(weight_1000_seeds_g/1000))%>%
  ggplot(aes(x=species,y=nb_seeds_needed))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  ylab("Nb seeds to get 0.5g (1mL macerate)")+
  geom_hline(aes(yintercept = 50),color="darkorange")+
  geom_hline(aes(yintercept = 20),color="darkorchid4")

data_quantity_seeds%>%
  group_by(common_name)%>%
  summarise(nb_leftovers_seeds=sum(nb_leftovers_seeds),
            leftover_weight_g=sum(leftover_weight_g))%>%
  left_join(select(data_ordered_seeds,common_name,family,genus,species,weight_1000_seeds_g))%>%
  arrange(desc(weight_1000_seeds_g))%>%
  mutate(species=forcats::fct_relevel(species,levels(data_ordered_seeds$species)))%>%
  ggplot(aes(x=species,y=nb_leftovers_seeds))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  geom_hline(aes(yintercept = 20),color="darkorchid4")+
  data_quantity_seeds%>%
  group_by(common_name)%>%
  summarise(nb_leftovers_seeds=sum(nb_leftovers_seeds),
            leftover_weight_g=sum(leftover_weight_g))%>%
  left_join(select(data_ordered_seeds,common_name,family,genus,species,weight_1000_seeds_g))%>%
  arrange(desc(weight_1000_seeds_g))%>%
  mutate(species=forcats::fct_relevel(species,levels(data_ordered_seeds$species)))%>%
  ggplot(aes(x=species,y=leftover_weight_g))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  geom_hline(aes(yintercept = 20),color="darkorchid4")

data_ordered_seeds%>%
  ggplot(aes(x=weight_1000_seeds_g))+
  geom_histogram(bins = 30)+
  data_ordered_seeds%>%
  ggplot(aes(x=weight_1000_seeds_g))+
  geom_histogram(bins = 30)+
  scale_x_log10()

# Insert in plant Megatree
test_tree <- rtrees::get_tree(sp_list = data_ordered_seeds,
                     taxon = "plant",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)
ggtree::ggtree(test_tree,
              layout="fan")+
  ggtree::geom_tiplab()

load("test2beignored/data/plant_phylogeny/plant_megatree.rda")
phytools::writeNexus(tree_plant_otl,"test2beignored/data/plant_phylogeny/nexus_plant.nxs")
# Leaf xploration

# install.packages("pliman")
library(pliman)
path_to_pics <- "../../Data/pics_greenhouse_and_leaves/leaves_4_area/cropped/"
leaves <- 
  image_import("P1190629_resultat.jpg",
               path = path_to_pics,
               plot = TRUE)
image_index(leaves,
            has_white_bg=T)
count_leaves <- analyze_objects(leaves,
                                fill_hull = T)


library(magick)
img <- image_read('../../Data/pics_greenhouse_and_leaves/leaves_4_area/cropped/P1190629_resultat.jpg')
img_bw <- img %>%
  image_convert(colorspace = "Gray") %>% 
  image_threshold(type = "black", threshold = "50%") %>%
  image_threshold(type = "white", threshold = "50%")
img_bw%>%
  image_threshold()



img_array <- as.integer(img_bw[[1]])
img_array[1,1,1]
img_array[1,1,2]
img_array[1,1,3]
