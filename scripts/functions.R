# FUNCTIONS #
# May 2020
# EKB

source('Scripts/find_millet_OTUs.R')

# READING IN FILES #============================================================

read_in_trnL_files <- function(path = "Data/SequencedData/Plants/RawData/PreppedFiles/trnL",
                               files = dir(path, pattern = "*.csv")){
  # this function reads all CSV files from set path and puts them in a list
  data <- files %>% 
    map(~ read_csv(file.path(path, .)))
  return(data)
}

read_in_ITS2_files <- function(path = "Data/SequencedData/Plants/RawData/PreppedFiles/ITS2",
                               files = dir(path, pattern = "*.csv")){
  # this function reads all CSV files from set path and puts them in a list
  data <- files %>% 
    map(~ read_csv(file.path(path, .)))
  return(data)
}

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

data_prep_multivariate <- function(list){
  
  # this function preps data for running NMDS
  # list argument is output from `filter_reads_data` fxn
  
  samples <- list[[1]]
  reads <- list[[3]]
  
  # convert to appropriate format for NMDS
  reads <- select(reads, SampleID, OTU, Rel_Reads)
  reads_spread <- pivot_wider(reads, names_from = OTU, values_from = Rel_Reads)
  reads_spread[is.na(reads_spread)] = 0
  reads_spread <- reads_spread %>% tibble::column_to_rownames("SampleID")
  
  sampleID <- intersect(reads$SampleID, samples$vial_barcode)
  groups <- samples %>% 
    filter(vial_barcode %in% sampleID)
  groups <- groups[match(sampleID, groups$vial_barcode),]
  
  return_list <- list(reads_spread, sampleID, groups)
  names(return_list) <- c("reads_spread", "sampleID", "groups")
  return(return_list)
  
}

run_metaMDS_til_converge <- function(df, prev_best, dist_metric, n) {
  repeat {
    # do something
    prev_best <- metaMDS(df, distance = dist_metric, 
                         trymax = 50, k = n, 
                         previous.best = prev_best)
    # exit if condition is met
    if (prev_best$converged == TRUE) break
  }
  return(prev_best)
}

binarize <- function(df){
  df[df > 0] <- 1
  return(df)
}

NMDS_plotting_prep <- function(data_list, dist_matrix) {
  
  # this function makes ellipses for plotting and gets perMANOVA values
  # data_list = output of `data_prep_multivariate` fxn
  
  # make dataframe with with MDS values and grouping variable
  groups <- data_list[[3]]
  NMDS <- data.frame(MDS1 = dist_matrix$points[,1], 
                     MDS2 = dist_matrix$points[,2], 
                     group = groups$group)
  
  # get mean point for each group
  NMDS.mean <- aggregate(NMDS[,1:2], list(group = groups$group), mean)
  
  # save results of ordiellipse() as an object
  plot(dist_matrix$points)
  ord <- ordiellipse(dist_matrix, 
                     groups$group, 
                     display = "sites", 
                     kind = "se", 
                     conf = 0.95, 
                     label = T)
  
  # make ellipses dataframe for ggplot2 
  df_ell <- data.frame()
  for(g in levels(NMDS$group)) {
    df_ell <-
      rbind(df_ell, cbind(as.data.frame(with(
        NMDS[NMDS$group == g,],
        vegan:::veganCovEllipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale)
      )),
      group = g))
  }
  
  # run perMANOVA
  group = as.matrix(groups$group)
  perMANOVA_output <- adonis(data_list[[1]] ~ group, permutations = 10000)
  
  # list of objects to return
  return_list <- list(NMDS, NMDS.mean, df_ell, perMANOVA_output)
  names(return_list) <- c("NMDS", "NMDS.mean", "df_ell", "perMANOVA_output")
  return(return_list)
  
}


# COMBINE DATA, REMOVE OUTLIERS, & PREP FOR PLOTTING #==========================

# TRNL # 

prep_460_allsp_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S010049")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S010049")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010049"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                    Factor = as.factor(data[[3]]$group),
                                    nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_466_allsp_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S013067")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S013067")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013067"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_460_PPonly_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  
  #remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S010049")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S010049")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010049"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_466_PPonly_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  
  #remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S013067")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S013067")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013067"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_454_PPonly_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()

  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_454_allsp_relabund <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_trnL(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}


# ITS2 #

prep_460_allsp_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S010044", "S010014", "S008810")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S010044", "S010014", "S008810")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010044", "S010014", "S008810"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_466_allsp_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S013043", "S013041")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S013043", "S013041")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013043", "S013041"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix,
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)

  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_460_PPonly_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S010044", "S010014", "S008810", "S010074")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S010044", "S010014", "S008810", "S010074")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010044", "S010014", "S008810", "S010074"),]
  
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_466_PPonly_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  
  # remove outliers
  data[[1]] <- 
    data[[1]][!(row.names(data[[1]]) %in% c("S013043", "S013041")),]
  data[[2]] <- 
    data[[2]][!data[[2]] %in% c("S013043", "S013041")]
  data[[3]] <- 
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013043", "S013041"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_454_PPonly_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}

prep_454_allsp_relabund_ITS2 <- function(samples, reads, totals, reads_min, period_code, rel_reads_min){
  
  data <- filter_reads_data_ITS2(samples, 
                                 reads, 
                                 totals, 
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate()
  data[[1]] <- binarize(data[[1]])
  
  if (rel_reads_min  <= 0.025){
    data[[1]] <- 
      data[[1]][!(row.names(data[[1]]) %in% c("S008824")),]
    data[[2]] <- 
      data[[2]][!data[[2]] %in% c("S008824")]
    data[[3]] <- 
      data[[3]][!(data[[3]]$vial_barcode) %in% c("S008824"),]
  } else {
    data[[1]] <- 
      data[[1]][!(row.names(data[[1]]) %in% c("S008824", "S008844", "S008847")),]
    data[[2]] <- 
      data[[2]][!data[[2]] %in% c("S008824", "S008844", "S008847")]
    data[[3]] <- 
      data[[3]][!(data[[3]]$vial_barcode) %in% c("S008824", "S008844", "S008847"),]
  }
  
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  dist_matrix <- metaMDSredist(dist_trnL)
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  
  pairwise_perMANOVA <- EcolUtils::adonis.pair(dist.mat = dist_matrix, 
                                               Factor = as.factor(data[[3]]$group),
                                               nper = 10000)
  
  # list of objects to return
  return_list <- list(df, pairwise_perMANOVA)
  names(return_list) <- c("df", "pairwise_perMANOVA")
  return(return_list)
  
}  

# FUNCTIONS FOR USING WEETUS #==================================================

data_prep_multivariate_WTU <- function(list){
  
  # this function preps data for running NMDS
  # list argument is output from `filter_reads_data` fxn
  
  samples <- list[[1]]
  reads <- list[[3]]
  
  # convert to appropriate format for NMDS
  reads <- select(reads, SampleID, WTU, Rel_Reads)
  reads_spread <- pivot_wider(reads, names_from = WTU, values_from = Rel_Reads)
  reads_spread[is.na(reads_spread)] = 0
  reads_spread <- reads_spread %>% tibble::column_to_rownames("SampleID")
  
  sampleID <- intersect(reads$SampleID, samples$vial_barcode)
  groups <- samples %>% 
    filter(vial_barcode %in% sampleID)
  groups <- groups[match(sampleID, groups$vial_barcode),]
  
  return_list <- list(reads_spread, sampleID, groups)
  names(return_list) <- c("reads_spread", "sampleID", "groups")
  return(return_list)
  
}

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

prep_466_allsp_relabund_WTU_genus <- function(samples, 
                                               reads, 
                                               totals, 
                                               OTU_WTU_key, 
                                               sum_taxa, 
                                               reads_min, 
                                               period_code, 
                                               rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                 reads, 
                                 totals, 
                                 OTU_WTU_key,
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_466_allsp_relabund_WTU_family <- function(samples, 
                                                reads, 
                                                totals, 
                                                OTU_WTU_key, 
                                                sum_taxa, 
                                                reads_min, 
                                                period_code, 
                                                rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S013017", "S013025")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S013017", "S013025")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013017", "S013025"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_466_PPonly_relabund_WTU_genus <- function(samples, 
                                          reads, 
                                          totals, 
                                          OTU_WTU_key, 
                                          sum_taxa, 
                                          reads_min, 
                                          period_code, 
                                          rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                 reads, 
                                 totals, 
                                 OTU_WTU_key,
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_466_PPonly_relabund_WTU_family <- function(samples, 
                                                reads, 
                                                totals, 
                                                OTU_WTU_key, 
                                                sum_taxa, 
                                                reads_min, 
                                                period_code, 
                                                rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  # remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S013017", "S013025")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S013017", "S013025")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013017", "S013025"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_460_allsp_relabund_WTU_genus <- function(samples, 
                                              reads, 
                                              totals, 
                                              OTU_WTU_key, 
                                              sum_taxa, 
                                              reads_min, 
                                              period_code, 
                                              rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_460_allsp_relabund_WTU_family <- function(samples, 
                                               reads, 
                                               totals, 
                                               OTU_WTU_key, 
                                               sum_taxa, 
                                               reads_min, 
                                               period_code, 
                                               rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S013017", "S013025")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S013017", "S013025")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S013017", "S013025"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_460_PPonly_relabund_WTU_genus <- function(samples, 
                                               reads, 
                                               totals, 
                                               OTU_WTU_key, 
                                               sum_taxa, 
                                               reads_min, 
                                               period_code, 
                                               rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  # remove outliers
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S010049")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S010049")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010049"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_460_PPonly_relabund_WTU_family <- function(samples, 
                                                reads, 
                                                totals, 
                                                OTU_WTU_key, 
                                                sum_taxa, 
                                                reads_min, 
                                                period_code, 
                                                rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       period_code = period_code, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_454_PPonly_relabund_WTU <- function(samples, 
                                          reads, 
                                          totals, 
                                          OTU_WTU_key, 
                                          sum_taxa,
                                          reads_min, 
                                          period_code, 
                                          rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                 reads, 
                                 totals, 
                                 OTU_WTU_key,
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_454_allsp_relabund_WTU <- function(samples, 
                                         reads, 
                                         totals, 
                                         OTU_WTU_key, 
                                         sum_taxa,
                                         reads_min, 
                                         period_code, 
                                         rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                 reads, 
                                 totals, 
                                 OTU_WTU_key,
                                 reads_min = reads_min, 
                                 period_code = period_code, 
                                 rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
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


prep_2017_allsp_relabund_WTU_ITS2 <- function(samples, 
                                              reads, 
                                              totals, 
                                              OTU_WTU_key, 
                                              sum_taxa, 
                                              reads_min, 
                                              yr, 
                                              rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       yr = yr, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  # remove outliers
  # group 1: "S008810", "S010014", "S013043"
  # group 2: group 1 + "S010063", "S010044", "S010012"
  data[[1]] <-
    data[[1]][!(row.names(data[[1]]) %in% c("S010049")),]
  data[[2]] <-
    data[[2]][!data[[2]] %in% c("S010049")]
  data[[3]] <-
    data[[3]][!(data[[3]]$vial_barcode) %in% c("S010049"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_2017_PPonly_relabund_WTU_ITS2 <- function(samples, 
                                               reads, 
                                               totals, 
                                               OTU_WTU_key, 
                                               sum_taxa,
                                               reads_min, 
                                               yr, 
                                               rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       yr = yr, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  # # remove outliers
  # data[[1]] <- 
  #   data[[1]][!(row.names(data[[1]]) %in% c("S010049", "S013067")),]
  # data[[2]] <- 
  #   data[[2]][!data[[2]] %in% c("S010049", "S013067")]
  # data[[3]] <- 
  #   data[[3]][!(data[[3]]$vial_barcode) %in% c("S010049", "S013067"),]
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_2016_PPonly_relabund_WTU_ITS2 <- function(samples, 
                                               reads, 
                                               totals, 
                                               OTU_WTU_key, 
                                               sum_taxa,
                                               reads_min, 
                                               yr, 
                                               rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       yr = yr, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}

prep_2016_allsp_relabund_WTU_ITS2 <- function(samples, 
                                              reads, 
                                              totals, 
                                              OTU_WTU_key,
                                              sum_taxa,
                                              reads_min, 
                                              yr, 
                                              rel_reads_min){
  
  data <- filter_reads_data_WeeTU_trnL(samples, 
                                       reads, 
                                       totals, 
                                       OTU_WTU_key,
                                       reads_min = reads_min, 
                                       yr = yr, 
                                       rel_reads_min = rel_reads_min) %>% 
    data_prep_multivariate_WTU()
  data[[1]] <- binarize(data[[1]])
  
  dist_trnL <- metaMDS(data[[1]], distance = "bray", trymax = 250, k = 3)
  plotting_data <- NMDS_plotting_prep(data, dist_trnL) 
  
  plotting_data[[1]]$df <- "NMDS"
  plotting_data[[2]]$df <- "NMDS.mean"
  plotting_data[[3]]$df <- "df_ell"
  plotting_data[[3]] <- plotting_data[[3]] %>% 
    rename("MDS1" = NMDS1, "MDS2" = NMDS2)
  df <- bind_rows(plotting_data[[1]], plotting_data[[2]], plotting_data[[3]])
  df$F.model <- plotting_data[[4]]$aov.tab$F.Model[1]
  df$pval <- plotting_data[[4]]$aov.tab$`Pr(>F)`[1]
  df$min_total <- reads_min
  df$min_rel_abund <- rel_reads_min
  df$sum_taxa <- sum_taxa
  
  return(df)
  
}
