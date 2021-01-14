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
krat_dat <- fecal %>% 
  filter(species == 'DO' | species == 'DM')

krat_trnl <- trnl %>% 
  filter(SampleID %in% krat_dat$vial_barcode)

krat_its <- its %>% 
  filter(SampleID %in% krat_dat$vial_barcode)

### remove millet OTUs

## trnL data
trnl_no_millet <- krat_trnl %>% 
  filter(!OTU %in% millet_OTUs$OTU)

# compare trnl before with millet removed
head(krat_trnl)
head(krat_trnl_no_millet)

## its data
krat_its_no_millet <- krat_its %>% 
  filter(!OTU %in% millet_OTUs_ITS2_no.hirt$OTU)

head(krat_its)
head(krat_its_no_millet)


