# Load packages and functions
library(dplyr)
library(magrittr)
library(ggplot2)
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
#   verbose = FALSE
# )

# Read file
data_ordered_seeds <- readxl::read_xlsx("test2beignored/data/monitoring_and_info.xlsx",sheet = 2) #list of ordered seeds
data_quantity_seeds <- readxl::read_xlsx("test2beignored/data/monitoring_and_info.xlsx",sheet = 5) #seeds sowing and stock tracking
data_phenotype_plants <- readxl::read_xlsx("test2beignored/data/phenotype.xlsx",sheet = 1) # phenotype info

# Format taxonomy
data_ordered_seeds%<>%
  select(!starts_with("...")) #remove empty columns generated when importing
data_ordered_seeds%<>% #rename some columns
  rename(source=`Achat - obtention des graines`,
         genotype_set=`Génotype/ Lot`,
         weight=Poids,
         n_seeds=`nombre de graines`,
         weight_1000_seeds=PMG)
data_ordered_seeds%<>% #create new taxonomy columns
  mutate(genus=stringr::str_extract(`Species for tree`,"^\\S*"),
         species=stringr::str_extract(`Species for tree`,"\\S*$"),
         family=Family)
data_ordered_seeds %<>%
  mutate(species=ifelse(nchar(species)<1,"sp",species))%>%
  mutate(species=paste0(genus,"_",species))
data_ordered_seeds%<>% #convert to numeric, will generate NA because some are text but its temporary
  mutate(n_seeds=as.numeric(n_seeds),
         weight_1000_seeds=as.numeric(weight_1000_seeds))%>%
  arrange(desc(weight_1000_seeds))%>%
  mutate(species=forcats::fct_inorder(species))

# Seed mass
data_ordered_seeds%>%
  ggplot(aes(x=species,y=weight_1000_seeds))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  ylab("1000 seeds weight (g)")+
  data_ordered_seeds%>%
  mutate(nb_seeds_needed=0.5/(weight_1000_seeds/1000))%>%
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
  rename(`Common name`=`espèce `,
         weight_1000_seeds=PMG)%>%
  left_join(select(data_ordered_seeds,`Common name`,species,family,genus))%>%
  arrange(desc(weight_1000_seeds))%>%
  mutate(species=forcats::fct_relevel(species,levels(data_ordered_seeds$species)))%>%
  ggplot(aes(x=species,y=`graines restantes`))+
  geom_point()+
  scale_y_log10()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.x = element_blank(),
        text = element_text(face="bold"))+
  geom_hline(aes(yintercept = 20),color="darkorchid4")

# Insert in plant Megatree


test_tree <- rtrees::get_tree(sp_list = data_ordered_seeds,
                     taxon = "plant",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)
ggtree::ggtree(test_tree,
              layout="fan")+
  ggtree::geom_tiplab()



# Leaf xploration

# install.packages("pliman")
library(pliman)
path_to_pics <- "../../Data/pics_greenhouse_and_leaves/"
leaves <- 
  image_import("P1190577.JPG",
               path = path_to_pics,
               plot = TRUE)
image_index(leaves)
count <- analyze_objects(leaves,
                         marker = "id",
                         fill_hull = T)
