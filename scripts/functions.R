# FUNCTIONS #
# May 2020
# EKB

source('Scripts/find_millet_OTUs.R')

# QUICK ADDITIONS #=============================================================

add_plot_type <- function(data){
  # add plot type to dataset -- control or KR exclosure
  data$plot_type <- NA
  for (i in 1:length(data$plot_type)) {
    if (data$plot[i] %in% c(4, 11, 14, 17)) {
      data$plot_type[i] = 'Control'
    } else {
      data$plot_type[i] = 'KR_Exclosure'
    }
  }
  return(data)
}

add_plotting_group <- function(data){
  # create grouping column based on species and plot
  data$group = NA                     
  for (i in 1:length(data$species)) {
    if (data$species[i] %in% c('DO', 'DM')) {
      data$group[i] = 'K-Rat'
    } else if (data$plot[i] %in% c(4, 11, 14, 17)) {
      data$group[i] = 'CP: Control'
    } else {
      data$group[i] = 'CP: KR Exclosure'
    }
  }
  return(data)
}

# SUMMARIZE DATA BY WEETU #=====================================================

summarize_trnL_by_WeeTU <- function(data, 
                                    col_quotes, 
                                    col_no_quotes){
  
  # make dataframe with all unique species
  all_species <- data %>% 
    subset(., !duplicated(WTU.species)) %>% 
    select(Kingdom:WTU.species)
  
  # select unique values to specified taxa level
  if (col_quotes %in% c("Species", "WTU.species")) {
    all_taxa <- all_species
  } else if (col_quotes %in% c("Genus", "WTU.genus")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species)
  } else if (col_quotes %in% c("Subfamily", "WTU.subfamily")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus)
  } else if (col_quotes %in% c("Family", "WTU.family")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Subfamily, -WTU.subfamily)
  } else if (col_quotes %in% c("Order", "WTU.order")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Subfamily, WTU.subfamily,
                                       -Family, -WTU.family)
  } else if (col_quotes %in% c("Clade2", "WTU.clade2")) {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom,
                                       Clade1, WTU.clade1,
                                       Clade2, WTU.clade2)
  } else if (col_quotes %in% c("Clade1", "WTU.clade1")) {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom, Clade1, WTU.clade1)
  } else  {
    all_taxa <- all_species %>% select(Kingdom, WTU.kingdom)
  }
  
  column <- enquo(col_no_quotes)
  
  sum <- data %>%
    group_by(SampleID, !! column) %>%
    summarise(Reads = sum(Reads))
  sum <- left_join(sum, all_taxa) %>% distinct()
  
  return(sum)
  
}


summarize_ITS2_by_WeeTU <- function(data, 
                                    col_quotes, 
                                    col_no_quotes){
  
  # make dataframe with all unique species
  all_species <- data %>% 
    subset(., !duplicated(WTU.species)) %>% 
    select(Domain:WTU.species)
  
  # select unique values to specified taxa level
  if (col_quotes %in% c("Species", "WTU.species")) {
    all_taxa <- all_species
  } else if (col_quotes %in% c("Genus", "WTU.genus")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species)
  } else if (col_quotes %in% c("Family", "WTU.family")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus)
  } else if (col_quotes %in% c("Order", "WTU.order")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       - Family, -WTU.family)
  } else if (col_quotes %in% c("Class", "WTU.class")) {
    all_taxa <- all_species %>% select(-Species, -WTU.species,
                                       -Genus, -WTU.genus,
                                       -Family, -WTU.family,
                                       -Order, -WTU.order)
  } else if (col_quotes %in% c("Clade1", "WTU.clade1")) {
    all_taxa <- all_species %>% select(Domain, WTU.domain, Clade1, WTU.clade1)
  } else  {
    all_taxa <- all_species %>% select(Domain, WTU.domain)
  }
  
  column <- enquo(col_no_quotes)
  
  sum <- data %>%
    group_by(SampleID, !! column) %>%
    summarise(Reads = sum(Reads))
  sum <- left_join(sum, all_taxa) %>% distinct()
  
  return(sum)
  
}

# PREP FOR NMDS #===============================================================

filter_reads_data_trnL <- function(samples,
                                   reads,
                                   totals,
                                   period_code = c(454, 460, 466),
                                   reads_min = 2000,
                                   rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), period %in% period_code) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !OTU %in% millet_OTUs$OTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}

filter_reads_data_ITS2 <- function(samples,
                                   reads,
                                   totals,
                                   period_code = c(454, 460, 466),
                                   reads_min = 2000,
                                   rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), period %in% period_code) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !OTU %in% millet_OTUs_ITS2_no.hirt$OTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}


# FUNCTIONS FOR USING WEETUS #==================================================

# TRNL # 

filter_reads_data_WeeTU_trnL <- function(samples,
                                         reads,
                                         totals,
                                         OTU_WTU_key,
                                         period_code = c(454, 460, 466),
                                         reads_min = 2000,
                                         rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), period %in% period_code) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # get millet WTUs/OTUs for taxa level
  millet_WTUs <- left_join(millet_OTUs, OTU_WTU_key)
  millet_WTUs <- millet_WTUs[, colnames(reads[2])] %>% 
    distinct() %>% 
    na.omit()
  colnames(millet_WTUs) <- c("WTU")
  
  # have WeeTU represent OTU
  reads <- reads[,c(1,3,2)]
  names(reads)[length(names(reads))] <- "WTU" 
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !WTU %in% millet_WTUs$WTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}


# ITS2 #

filter_reads_data_WeeTU_ITS2 <- function(samples,
                                         reads,
                                         totals,
                                         OTU_WTU_key,
                                         period_code = c(454, 460, 466),
                                         reads_min = 2000,
                                         rel_reads_min = 0.001){
  
  # add plot type to fecal collection data
  # add group for plotting
  # and remove samples that were part of the trap/bait test
  samples <- add_plot_type(samples) %>% 
    add_plotting_group() %>% 
    filter(is.na(notes), period %in% period_code) 
  
  # select only fecal samples
  fecal_id <- samples$vial_barcode
  
  # get millet WTUs/OTUs for taxa level
  millet_WTUs <- left_join(millet_OTUs_ITS2_no.hirt, OTU_WTU_key)
  millet_WTUs <- millet_WTUs[, colnames(reads[2])] %>% 
    distinct() %>% 
    na.omit()
  colnames(millet_WTUs) <- c("WTU")
  
  # have WeeTU represent OTU
  reads <- reads[,c(1,3,2)]
  names(reads)[length(names(reads))] <- "WTU" 
  
  # add totals to reads df
  # select only fecal samples and remove millet OTUs
  # and make relative reads column
  reads <- full_join(reads, totals)
  reads <- reads %>% 
    filter(SampleID %in% fecal_id, !WTU %in% millet_WTUs$WTU) %>% 
    mutate(Rel_Reads = Reads/Total_Reads)
  
  # filter data by minimum total reads and/or minimum relative reads
  reads <- reads %>% 
    filter(Total_Reads >= reads_min, Rel_Reads >= rel_reads_min)
  
  return_list <- list(samples, fecal_id, reads)
  names(return_list) <- c("samples", "fecal_id", "reads")
  return(return_list)
  
}

