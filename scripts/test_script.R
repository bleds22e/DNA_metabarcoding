# TEST SCRIPT #

# LIBRARIES #
library(tidyverse)
source("scripts/functions.R")

# DATA #

# DNA data
trnl <- read_csv("data/sequenced_data/trnL_reads_WeeTU.csv")
its <- read_csv("data/sequenced_data/ITS2_reads_WeeTU.csv")

# collection data 
fecal <- read_csv("data/collection_data/fecal_sample_collection.csv")

### select SampleID from k-rats only
# also need to select just data from traps and without any added notes
krat_dat <- fecal %>% 
  filter(species == 'DO' | species == 'DM', 
         sample_type == 'trap', 
         is.na(notes)) 
write_csv(krat_dat, "data/collection_data/krat_fecal_data.csv")

krat_trnl <- trnl %>% 
  filter(SampleID %in% krat_dat$vial_barcode)
write_csv(krat_trnl, "data/sequenced_data/krat_fecal_trnl.csv")

krat_its <- its %>% 
  filter(SampleID %in% krat_dat$vial_barcode)
write_csv(krat_its, "data/sequenced_data/krat_fecal_its.csv")

### remove millet OTUs

## trnL data
krat_trnl_no_millet <- krat_trnl %>% 
  filter(!OTU %in% millet_OTUs$OTU)
write_csv(krat_trnl_no_millet, "data/sequenced_data/krat_fecal_trnl_no_millet.csv")

# compare trnl before with millet removed
head(krat_trnl)
head(krat_trnl_no_millet)

## its data
krat_its_no_millet <- krat_its %>% 
  filter(!OTU %in% millet_OTUs_ITS2_no.hirt$OTU)
write_csv(krat_its_no_millet, "data/sequenced_data/krat_fecal_its_no_millet.csv")

head(krat_its)
head(krat_its_no_millet)

###
# if you want to look at data from data from one specific period, you can use
# the krat_dat dataframe to filter which samples come from which period
