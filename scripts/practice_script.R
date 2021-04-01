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
         Domain == "Eukaryota")

# find the average number of reads per plant family for each krat species
krat_family_reads <- krat_sp_its %>%          # use the krat_sp_its df
  group_by(rodent_sp, Family, period) %>%     # group the data by these 3 columns
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
# find the average number of reads per plant family for each krat species
krat_Genus_reads <- krat_sp_its %>%          # use the krat_sp_its df
  group_by(rodent_sp, Genus, period) %>%     # group the data by these 3 columns
  summarize(avg_genus_reads = mean(Reads))   # make a new column called avg_family_reads
# use the package ggplot2 to make some plots
ggplot(data = krat_Genus_reads, aes(x = Genus, y = avg_genus_reads)) +
  # tell ggplot what df to use and which columns will be the values on the x and y axes
  geom_col() + # add a geom (like a shape) to plot--this one is saying make a plot with columns
  facet_wrap(period ~ rodent_sp,     # split the plot into multiple plots, making a mini plot for 
             nrow = 3, ncol = 2) +   # each combination of period and krat. Make it have 3 rows and 2 columns
  theme_bw() + # make it prettier
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

# Other things that you might want to investigate or think about:
# - remove any OTUs that have fewer than 10 or 50 or 100 reads to focus on the majority of the diet
# - use presence/absense and see which plants DMs and DOs have in common
# - compare what the krats are eating to what plants we find on the control plots
# - use a diet selectivity metric to see how what plants DMs and DOs select compare to 
#   what plants are available on each plot and if one species is more selective than the other.
#   One metric that I know of is Chesson's index, but I've never done this before.