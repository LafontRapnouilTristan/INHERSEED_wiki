# RTRY testing
# need manual requests... :/
install.packages("rtry")
library(rtry)

#get list of TRY species
download.file('https://try-db.org/dnld/TryAccSpecies.txt', destfile = "test2beignored/data/TryAccSpecies.txt", method = "wget", extra = "-r -p --random-wait")

try_sp <- data_ordered_seeds %>%
  mutate(scientific_cleaned = gsub("_"," ",species))%>%
  select(scientific_cleaned) %>%
  mutate(scientific_cleaned=ifelse(scientific_cleaned=="Legousia veneris","Legousia speculum-veneris",scientific_cleaned))%>%
  na.omit() %>%
  unique() %>%
  left_join(readr::read_tsv("test2beignored/data/TryAccSpecies.txt"),
            by = c("scientific_cleaned" = "AccSpeciesName")
  )


try_request <- try_sp%>%
  summarise(species_list=paste0(AccSpeciesID,collapse = ", "))
try_request%>%
  readr::write_tsv("test2beignored/data/TryReqSpecies.txt")

cantarel_TRY <- readxl::read_xlsx("../../Data/list_traits_TRY/Traits repro TRY database.xlsx")
try_request_traits <- cantarel_TRY%>%
  filter(tokeep=="y")%>%
  select(TraitID)%>%
  summarise(traits_list=paste0(TraitID,collapse = ", "))
try_request_traits%>%
  readr::write_tsv("test2beignored/data/TryReqTraits.txt")
# 41259 - ID of the first request 



owo <- readr::read_tsv("../../Data/list_traits_TRY/41259_1905202519280_TRY_relaease/41259.txt")
owo %<>%
  select(AccSpeciesName, TraitName, StdValue) %>%
  rename(
    species = AccSpeciesName,
    trait = TraitName,
    value = StdValue
  ) %>%
  group_by(species, trait) %>%
  summarise(value = mean(as.numeric(value), na.rm = TRUE)) %>%
  tidyr::pivot_wider(values_from = "value", names_from = "trait")%>%
  mutate(across(everything(),~ifelse(.x=="NaN",NA,.x)))



owo%>%
  reshape2::melt(id.vars=c("species"))%>%
  group_by(variable)%>%
  mutate(is_na_var=ifelse(is.na(value),0,1))%>%
  summarise(pct_species_avail=(sum(is_na_var)/82)*100)%>%
  ggplot(aes(x=pct_species_avail,y=variable))+
  geom_point()+
  theme_minimal()



owo%>%
  reshape2::melt(id.vars=c("species"))%>%
  group_by(variable)%>%
  mutate(is_na_var=ifelse(is.na(value),0,1))%>%
  summarise(pct_species_avail=(sum(is_na_var)/82)*100)%>%
  filter(pct_species_avail>0)%>%
  readr::write_tsv("../../Data/list_traits_TRY/traits_with_some_values.txt")
 