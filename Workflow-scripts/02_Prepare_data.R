
#
# 01
# Data preparation
#



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Load data image ----
load("Data/RData/00.RData")



# ---- Define sample information for final datasets ---- 

## Define groups
get_observations_data(which = "observations", 
                      output.type = "list_",
                      dataset = "proteinGroups") %>%
  do_mutate("groups", factor(.strsplit_keep(observations, "_", 2:3))) %>%
  do_fun(save_observations_data,
         column = "groups",
         name = "groups",
         dataset = "proteinGroups")



# ---- Clean variables ----

## Identify contminant and other proteins to be removed
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
save.image("Data/RData/01.RData")
