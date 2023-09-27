
#
# 02
# Data preparation
#



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Load data image ----
load("Data/RData/01.RData")



# ---- Define sample information for raw datasets ---- 


## Define technical replicates in data
get_observations_data(which = "observations") %>%
  do_mutate("replicates", factor(.strsplit_keep(observations, "_", 1:3))) %>%
  do_fun(save_observations_data,
         column = "replicates",
         name = "replicates", 
         dataset = "proteinGroups")




## Combine technical replicates based on Unique.peptides data 
## and make new dataset
get_data_frame(which = "Unique.peptides",
               variables = ,
               observations = ) %>%
  add_observations_data(which = "replicates") %>%
  do_row_summary(FUN = mean, by = "replicates") %>%
  do_fun(data_frame2new_dataset,
         "Mean.peptides",
         output.name = "New_dataset") %>%
  do_fun(save_dataset, 
         name = "proteinGroups_c")
  
  ## Combine LFQ intensity data
  get_data_frame(which = "LFQ.intensity",
                 variables = ,
                 observations = , 
                 dataset = "proteinGroups") %>%
  add_observations_data(which = "replicates") %>%
  do_row_summary(FUN = mean_or_0, by = "replicates") %>%
  do_fun(save_data_frame, 
         name = "Mean.LFQ.intensity", 
         dataset = "proteinGroups_c")


## Dataset Control
  get_observations_data(which = "observations", 
                        dataset = "proteinGroups_c") %>%
    do_mutate("genotype", factor(.strsplit_keep_first(observations, "_"))) %>% 
    do_fun(save_observations_data,
           column = "genotype",
           name = "genotype",
           dataset = "proteinGroups_c")



# ---- Transfer variables data from raw to final datasets ----
get_variables_data(which = c("Gene.names",
                             "Protein.names",
                             "Only.identified.by.site",
                             "Potential.contaminant",
                             "Reverse"),
                   dataset = "proteinGroups") %>%
  do_fun(save_variables_data, 
         column = c("Gene.names",
                    "Protein.names",
                    "Only.identified.by.site",
                    "Potential.contaminant",
                    "Reverse"),
         name = c("Gene.names",
                  "Protein.names",
                  "Only.identified.by.site",
                  "Potential.contaminant",
                  "Reverse"),
         dataset = "proteinGroups_c")


# ---- Clean variables ----

## Dataset Control
get_variables_data(which = c("Potential.contaminant",
                             "Only.identified.by.site",
                             "Reverse"),
                   dataset = "proteinGroups_c") %>%
    do_expr(FUN = is.na, modify = c("Potential.contaminant",
                                    "Only.identified.by.site",
                                    "Reverse")) %>%
  do_column_summary(FUN = function(x) mean(x) == 1,
                    name = "clean",
                    output.name = "variables_data_expr_columns") %>%
  do_fun(save_variables_data,
         column = "clean",
         dataset = "proteinGroups_c")



# ---- Save data image ----
save.image("Data/RData/02.RData")
