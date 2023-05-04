
#
# 00
# Import of raw data, trait information, and fasta files
#



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Create data structure ----
Info <- list(Imports = list())
Datasets <- list()
Analysis <- list()



# ---- Import data ----
import_files(files = c("Data/proteinGroups.txt"))



# ---- Dataset combined MQ search ----

# Dataset1
Datasets[["Dataset1"]] <-
  import2new_dataset(raw.data = Info$Imports[["proteinGroups"]],
                     variable.identifiers = 
                       .strsplit_keep_first(`Protein.IDs`, split = ";"),
                     variables.data = c("Gene.names",
                                        "Protein.names",
                                        "Only.identified.by.site",
                                        "Potential.contaminant",
                                        "Reverse"),
                     observations = c("A",
                                      "B",
                                      "C"),
                     data.frames = c("Unique.peptides",
                                     "LFQ.intensity"))



# ---- Save data image ----
save.image("Data/RData/00.RData")
