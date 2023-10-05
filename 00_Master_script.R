
#
# 00 
# Master script
#



# ---- Setup folders ---- 
dir.create("Data")
dir.create("Data/raw")
dir.create("Data/RData")
dir.create("Scripts")
dir.create("Output")



# ---- QC of raw files ---- 
franksQC("Data/PTM/")



# ---- 01_Import_data ---- 

## Collect file names
list.files("Data", recursive = T, pattern = "\\.") %>% 
  paste0("Data/", .) %>% 
  .cat_character()

## Identify samples names 
.identify_observations(Info[["Imports"]][[6]]) %>% 
  .cat_character()





# ---- Observation names ---- 
sample_names <- Info[["Imports"]][["Aligned_feature_list_filtered_min2"]] %>% 
  select(matches("area")) %>% 
  names() %>% 
  `[`(-1) %>% 
  gsub("datafile:", "", .) %>% 
  substring(1, regexpr(":", .) - 14)

.cat_character(sample_names)

rm(sample_names)



# ---- 02_Prepare_data ---- 

get_observations()





