
#
# 00 
# Master script
#



# ---- Check packages ---- 
packages <- c("tidyverse", 
              "BiocManager",
              "eulerr", 
              "ggrepel", 
              "ggvenn", 
              "plotly", 
              "Peptides")

(packages_installed <- 
    setNames(object = packages %in% rownames(installed.packages()), 
             nm = packages))

# Install packages from CRAN
for (i in 1:length(packages_installed)) {
  
  if (!packages_installed[[i]])
    install.packages(names(packages_installed)[[i]])
  
}



# Bioconductor packages
packages_bc <- c("clusterProfiler", 
                 "ComplexHeatmap", 
                 "org.Hs.eg.db", 
                 "org.Mm.eg.db", 
                 "GO.db")

(packages_bc_installed <- 
    setNames(object = packages_bc %in% rownames(installed.packages()), 
             nm = packages_bc))

# Install packages from CRAN
for (i in 1:length(packages_bc_installed)) {
  
  if (!packages_bc_installed[[i]])
    BiocManager::install(names(packages_bc_installed)[[i]])
  
}






# ---- Setup folders ---- 
dir.create("Data")
dir.create("Data/RData")
dir.create("Scripts")
dir.create("Output")





# ---- Download template files ---- 

## Check available templates
browseURL("https://github.com/nicohuttmann/pOmics2-templates/tree/master")

## List files to download
files <- c("01_Import_data.R",
           "02_Prepare_data.R")

## Set base url for workflow scripts
base_url <- "https://raw.githubusercontent.com/nicohuttmann/pOmics2-templates/master/Workflow-scripts/"

# Download files to folder Scripts
for (i in files) {
  download.file(url = paste0(base_url, i),
                destfile = paste("Scripts/", i),
                quiet = TRUE)
}

rm(list = c("files", 
            "base_url", 
            "i"))



# ---- Scripts ---- 
source("Scripts/01_Import_data.R")
source("Scripts/02_Prepare_data.R")
source("Scripts/03_Qualitative_analysis.R")



# ---- 01_Import_data ---- 


## Collect file names
list.files("Data", recursive = T, pattern = "\\.txt") %>% 
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










#' Guesses observations names from raw_dataset
#'
#' @param x data frame
#' @param pattern separation pattern of column names
#' @param sample.as.suffix
#' @param to.include observations that must be included; helps to identify
#' correct observations
#'
#' @return
#' @export
#'
#'
.identify_observations <- function(x,
                                   pattern = ".",
                                   sample.as.suffix = T,
                                   to.include = NULL) {
  
  if (sample.as.suffix) {
    fix <- substring(text = colnames(x), 1, .str_locate_last(colnames(x), pattern = pattern) - 1)
  } else {
    fix <- substring(text = colnames(x), .str_locate_last(colnames(x), pattern = pattern) + 1)
  }
  
  
  tab <- sort(table(fix), decreasing = TRUE)
  
  tab <- tab[nchar(names(tab)) >= 3]
  
  # Important: Exclusion list for potential sample names
  tab <- tab[!names(tab) %in% c("Count", "IDs", "acid", "window", "position",
                                "names", "[%]", "Peptide.counts")]
  
  # Consider predefined observations
  if(is.null(to.include)) {
    # Choose observations
    observations <- colnames(x)[grepl(names(tab)[1], colnames(x))] %>%
      substring(first = nchar(names(tab)[1]) + 2)
  } else {
    observations <- colnames(x)[grepl(names(tab)[1], colnames(x))] %>%
      substring(first = nchar(names(tab)[1]) + 2)
    while (!all(to.include) %in% observations) {
      if (length(tab) == 0) {
        stop(paste0(
          "Given observations ",
          paste(to.include, collapse = ", "),
          " could not all be identified in the data frame. Check your input or the algorithm."
        ))
      }
      tab <- tab[-1]
      observations <- colnames(x)[grepl(names(tab)[1], colnames(x))] %>%
        substring(first = nchar(names(tab)[1]) + 2)
    }
  }
  
  return(observations)
  
}


#' Locates last position of pattern in given string
#'
#' @param string vector of strings
#' @param pattern pattern to match
#'
#' @return
#' @export
#'
#'
.str_locate_last <- function(string, pattern) {
  
  pattern <- .str_rev(pattern)
  
  pattern <- gsub(pattern = "\\.", replacement = "\\\\.", x = pattern)
  
  
  pos <- regexpr(pattern = pattern, text = .str_rev(string))
  
  for (i in seq_along(pos)) {
    
    if (pos[i] != -1) pos[i] <- nchar(string)[i] - pos[i] + 1
    
  }
  
  
  return(c(pos))
  
}


#' Returns given strings in reverse order
#'
#' @param strings vector of strings
#'
#' @return
#' @export
#'
#'
.str_rev <- function(strings) {
  
  return(
    sapply(X = strings, FUN = function(x)
    {
      paste(rev(.strsplit(x, "")), collapse = "")
    }, USE.NAMES = F)
  )
  
}




