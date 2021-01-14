# Find millet OTUs to remove
# EKB
# 27 May 2020

# LIBRARIES #
library(tidyverse)

### trnL ###

reads_WTU <- read_csv("data/sequenced_data/trnL_reads_WeeTU.csv")

# create OTU_WTU key for trnL
all_OTU_WTU <- reads_WTU %>% 
  select(OTU, Kingdom:WTU.species) %>% 
  distinct()
#write_csv(all_OTU_WTU, "Data/SequencedData/Plants/ProcessedData/trnL_OTU_WTU_key.csv")

# find millet OTUs to remove
all_species <- reads_WTU %>% 
  subset(., !duplicated(WTU.species)) %>% 
  select(Kingdom:WTU.species)

millet_WTU.genus <- all_species %>% 
  filter(Family == "Poaceae", 
         Subfamily == "Panicoideae" | is.na(Subfamily),
         Genus == "Panicum" | is.na(Genus)) %>% 
  select(WTU.genus) %>% 
  distinct()

millet_OTUs <- reads_WTU %>% 
  filter(WTU.genus %in% millet_WTU.genus$WTU.genus) %>% 
  select(OTU) %>% 
  distinct()

### ITS2 ###

reads_WTU_ITS2 <- read_csv("data/sequenced_data/ITS2_reads_WeeTU.csv")

# create OTU_WTU key for trnL
all_OTU_WTU_ITS2 <- reads_WTU_ITS2 %>% 
  select(OTU, Domain:WTU.species) %>% 
  distinct()
#write_csv(all_OTU_WTU_ITS2, "Data/SequencedData/Plants/ProcessedData/ITS2_OTU_WTU_key.csv")

# find millet OTUs to remove
all_species_ITS2 <- reads_WTU_ITS2 %>% 
  subset(., !duplicated(WTU.species)) %>% 
  select(Domain:WTU.species)

# WTU.species 119:125 are Panicum
#  - WTU.species.122 is Panicum miliaceum
#  - WTU.species.121 is Panicum hirticaule
millet_WTU.sp_ITS2 <- all_species_ITS2 %>% 
  filter(Family == "Poaceae", 
         Genus == "Panicum" | is.na(Genus), 
         WTU.species != 121) %>% 
  select(WTU.species) %>% 
  distinct()

millet_OTUs_ITS2_no.hirt <- reads_WTU_ITS2 %>% 
  filter(WTU.species %in% millet_WTU.sp_ITS2$WTU.species) %>% 
  select(OTU) %>% 
  distinct()

millet_OTUs_ITS2_pani.mili.only <- reads_WTU_ITS2 %>% 
  filter(WTU.species == 122) %>% 
  select(OTU) %>% 
  distinct()

# no.hirt and pani.mili.only are almost identical (6095 OTUs vs. 5984 OTUs)
# probably best to be conservative, so will go with no.hirt

rm(reads_WTU, reads_WTU_ITS2, all_OTU_WTU, all_OTU_WTU_ITS2, 
   all_species, all_species_ITS2, millet_WTU.genus, millet_WTU.sp_ITS2,
   millet_OTUs_ITS2_pani.mili.only)
