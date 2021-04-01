# Practice Script #

library(tidyverse)

# load the data
# make sure the file paths match up
krat_its <- read_csv("data/sequenced_data/krat_fecal_its_no_millet.csv")
krat_fecal <- read_csv("data/collection_data/krat_fecal_data.csv")

# use a left_join function to add data from the krat_fecal df to the krat_its df
# definitely check out the other dplyr joins (google dplyr cheatsheet), as they
# can also be very helpful
krat_sp_its <- left_join(krat_its, 
                         select(krat_fecal, vial_barcode, species, period), 
                         # the line above makes a df with 3 columns from krat_fecal
                         # to merge with the krat_its df
                         by = c("SampleID" = "vial_barcode")) # this says match up each row
                              # of the two dfs we are joining based on the following columns
                              # because the columns have different names in each df, we have 
                              # to use the c() and = to say that they are equivalent columns

# rename the species column that we added above to be "rodent_sp"
krat_sp_its <- krat_sp_its %>% 
  rename("rodent_sp" = "species") %>% 
  filter(Reads >= 5,
         Domain == "Eukaryota", 
         !is.na(Genus))

all_samples <- krat_sp_its %>%        
  expand(SampleID, Genus)

# find the average number of reads per plant family for each krat species
krat_family_reads <- krat_sp_its %>% 
  full_join(all_samples) %>%                  # use the krat_sp_its df 
  group_by(rodent_sp, Family, period) %>% 
  replace_na(list(0)) %>%                                  # group the data by these 3 columns
  summarize(avg_family_reads = mean(Reads))   # make a new column called avg_family_reads
                                              # by taking the mean of the Reads column for 
                                              # for each group that we made above 
                                              # (i.e., Poaceae from DM caught in period 454)

# use the package ggplot2 to make some plots
ggplot(data = krat_family_reads, aes(x = Family, y = avg_family_reads)) +
      # tell ggplot what df to use and which columns will be the values on the x and y axes
  geom_col() + # add a geom (like a shape) to plot--this one is saying make a plot with columns
  facet_wrap(period ~ rodent_sp,     # split the plot into multiple plots, making a mini plot for 
             nrow = 3, ncol = 2) +   # each combination of period and krat. Make it have 3 rows and 2 columns
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) # rotate the family names on the x-axis


### What about with Genus?
# find the average number of reads per plant genus for each krat species
krat_genus_reads <- krat_sp_its %>% 
  full_join(all_samples) %>%                  
  group_by(rodent_sp, Genus, period) %>% 
  replace_na(list(0)) %>%  
  mutate(Genus=replace(Genus,Genus=="Hoffmannseggia","Hoffmanseggia")) %>%
  summarize(avg_genus_reads = mean(Reads)) %>%
  drop_na() %>%
  mutate(type = "Rodents")

# create color palate for ggplot to use
pal <- c(
  "DM" = "salmon",
  "DO" = "#4E84C4", 
  "abundance" = "#52854C"
  )

# Genus-level plots
ggplot(data = krat_genus_reads, aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  # tell ggplot what df to use and which columns will be the values on the x and y axes
  geom_col(width = .5, position = position_dodge(1)) + # add a geom (like a shape) to plot--this one is saying make a plot with columns
  facet_wrap(~period, nrow=3) +   
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))

# Bring in portal plant and rodent data
tables <- portalr::load_plant_data()
rodents <- portalr::load_rodent_data()

plants <- portalr::plant_abundance(level="treatment") %>%
  filter(treatment=="control") %>%
  left_join(tables$species_table) %>%
  filter(year %in% c(2016, 2017)) %>%
  filter(!(year == 2016 & season == "winter")) %>%
  mutate(period = 454, Genus=genus, type="Plants") %>%
  mutate(period = ifelse(year==2017 & season=="winter", 460, period)) %>%
  mutate(period = ifelse(year==2017 & season=="summer", 466, period)) %>%
  mutate(rodent_sp = "abundance") %>%
  select(rodent_sp, Genus, period, abundance, type) %>%
  group_by(rodent_sp, period, Genus, type) %>%
  summarise(avg_genus_reads = sum(abundance))

# Combine dna and plant tables, add season variable for plotting

rodents_plants <- bind_rows(krat_genus_reads, plants) %>%
  mutate(season = "Summer 2016") %>%
  mutate(season = replace(season, period==460, "Winter 2017"),
         season = replace(season, period==466, "Summer 2017"))

# Plot all data

ggplot(data = rodents_plants, aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  geom_col(width = .5, position = position_dodge(1)) +
  facet_wrap(~season, nrow=3) +   
  scale_y_log10(name = "count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))


# Focus on some highlights
# Asynchrony
rodents_plants %>%
  filter(Genus %in% c("Astragalus", "Amaranthus", "Bouteloua")) %>%
ggplot(aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  geom_col(position = "dodge") +
  facet_wrap(~season, nrow=3) +   
  scale_y_log10(name = "count") +
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(size = 10)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))

# Legume family
rodents_plants %>%
  filter(Genus %in% c("Astragalus","Vachellia","Dalea","Hoffmanseggia","Prosopis","Lupinus")) %>%
  ggplot(aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  geom_col(position = "dodge") + 
  facet_wrap(~season, nrow=3) +   
  scale_y_log10(name = "count") +
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(size = 9)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))

# Grasses
rodents_plants %>%
  filter(Genus %in% c("Aristida","Bouteloua","Eragrostis","Muhlenbergia","Munroa","Setaria",     
                      "Eriochloa","Panicum","Urochloa","Hordeum","Enneapogon","Digitaria","Hilaria",     
                      "Yakirra","Coleataenia","Chondrosum","Bromus","Sporobolus","Tragus","Digitaria",
                      "Erioneuron","Vulpia","Schismus","Lycurus","Chloris","Enneapogon","Buchloe")) %>%
  ggplot(aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  geom_col(position = "dodge") + 
  facet_wrap(~season, nrow=3) +   
  scale_y_log10(name = "count") +
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 9)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))

# DO Favorites
rodents_plants %>%
  filter(Genus %in% c("Plantago","Xanthisma","Nuttallanthus","Linanthus", "Kallestromia", "Eschscholzia",
                      "Euthamia")) %>%
  ggplot(aes(x = Genus, y = avg_genus_reads, fill = rodent_sp)) +
  geom_col(position = "dodge") + 
  facet_wrap(~season, nrow=3) +   
  scale_y_log10(name = "count") +
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 9)) +
  scale_fill_manual(name = "" , values = pal, limits = names(pal))
  

# Other things that you might want to investigate or think about:
# - use presence/absense and see which plants DMs and DOs have in common
# - compare what the krats are eating to what plants we find on the control plots
# - use a diet selectivity metric to see how what plants DMs and DOs select compare to 
#   what plants are available on each plot and if one species is more selective than the other.
#   One metric that I know of is Chesson's index, but I've never done this before.