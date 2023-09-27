
#
# 01
# Data preparation
#



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Load data image ----
load("Data/RData/01.RData")



# ---- Define sample information for raw datasets ---- 

## Combine dataset All_raw (Unique.peptides, LFQ.intensity)
Analysis[["Process_All_dataset"]] <- get_observations_data(which = "observations", 
                                                           output.type = "list_",
                                                           dataset = "proteinGroups_raw") %>%
  do_mutate("replicates", factor(.strsplit_keep(observations, "_", 1:3))) %>%
  do_fun(save_observations_data,
         column = "replicates",
         name = "replicates",
         dataset = "proteinGroups_raw") %>% 
  get_data_frame(which = "Unique.peptides",
                 variables = ,
                 observations = ,
                 output.type = "list_",
                 output.name = ,
                 dataset = "proteinGroups_raw") %>%
  add_observations_data(which = "replicates",
                        dataset = "proteinGroups_raw") %>%
  do_row_summary(FUN = mean, by = "replicates") %>%
  do_fun(data_frame2new_dataset,
         "Mean.peptides",
         output.name = "New_dataset") %>%
  do_fun(save_dataset, name = "proteinGroups") %>% 
  get_data_frame(which = "LFQ.intensity",
                 variables = ,
                 observations = ,
                 output.type = "list_",
                 output.name = ,
                 dataset = "proteinGroups_raw") %>%
  add_observations_data(which = "replicates",
                        dataset = "proteinGroups_raw") %>%
  do_row_summary(FUN = pOmics::mean_or_0, by = "replicates") %>%
  do_fun(save_data_frame, name = "Mean.LFQ.intensity", dataset = "proteinGroups")



# ---- Define sample information for final datasets ---- 

## Dataset All
Analysis[["Definitions_All"]] <- get_observations_data(which = "observations", 
                                                       output.type = "list_",
                                                       dataset = "proteinGroups") %>%
  do_mutate("dataset", factor(ifelse(observations %in% 
                                       get_observations(dataset = "Control"), 
                                     "Control", "Diet"))) %>% 
  do_fun(save_observations_data,
         column = "dataset",
         name = "dataset",
         dataset = "proteinGroups") %>%
  do_mutate("replicates", factor(.strsplit_keep(observations, "_", 1:3))) %>%
  do_fun(save_observations_data,
         column = "replicates",
         name = "replicates",
         dataset = "proteinGroups") %>% 
  do_mutate("genotype", factor(.strsplit_keep_first(observations, "_"))) %>% 
  do_fun(save_observations_data,
         column = "genotype",
         name = "genotype",
         dataset = "proteinGroups") %>% 
  do_mutate("diet", factor(.strsplit_keep(observations, "_", 2))) %>% 
  do_fun(save_observations_data,
         column = "diet",
         name = "diet",
         dataset = "proteinGroups") %>% 
  do_mutate("group", factor(.strsplit_keep(observations, "_", 1:2))) %>% 
  do_fun(save_observations_data,
         column = "group",
         name = "group",
         dataset = "proteinGroups")



# ---- Transfer variables data from raw to final datasets ----

## Dataset All
get_variables_data(which = c("Gene.names",
                             "Protein.names",
                             "Only.identified.by.site",
                             "Potential.contaminant",
                             "Reverse"),
                   dataset = "proteinGroups_raw") %>%
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
         dataset = "proteinGroups")



# ---- Clean variables ----

## Dataset All
get_variables_data(which = c("Potential.contaminant",
                             "Only.identified.by.site",
                             "Reverse"),
                   dataset = "proteinGroups") %>%
  do_expr(is.na(x), modify = c("Potential.contaminant",
                               "Only.identified.by.site",
                               "Reverse")) %>%
  do_column_summary(FUN = function(x) mean(x) == 1,
                    name = "clean",
                    output.name = "variables_data_expr_columns") %>%
  do_fun(save_variables_data,
         column = "clean",
         dataset = "proteinGroups")



# ---- Save data image ----
save.image("Data/RData/02.RData")
