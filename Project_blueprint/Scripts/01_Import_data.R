
#
# 01
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
import_files(files = "Data/proteinGroups.txt")



# ---- Tidy up data sets ----
Datasets[["proteinGroups"]] <-
  import2new_dataset(raw.data = Info$Imports[["proteinGroups"]],
                     variable.identifiers = .strsplit_keep_first(`Protein.IDs`, split = ";"),
                     variables.data = c("Gene.names",
                                        "Protein.names",
                                        "Only.identified.by.site",
                                        "Potential.contaminant",
                                        "Reverse"),
                     observations = c("IRP1_Ct_03_A",
                                      "IRP1_Ct_03_B",
                                      "IRP1_Ct_04_A",
                                      "IRP1_Ct_04_B",
                                      "IRP1_Ct_05_A",
                                      "IRP1_Ct_05_B",
                                      "IRP1_Ct_06_A",
                                      "IRP1_Ct_06_B",
                                      "IRP1_Ct_07_A",
                                      "IRP1_Ct_07_B",
                                      "WT_Ct_2279_A",
                                      "WT_Ct_2279_B",
                                      "WT_Ct_2280_A",
                                      "WT_Ct_2280_B",
                                      "WT_Ct_2281_A",
                                      "WT_Ct_2281_B",
                                      "WT_Ct_2282_A",
                                      "WT_Ct_2282_B",
                                      "WT_Ct_2283_A",
                                      "WT_Ct_2283_B"),
                     data.frames = c("Unique.peptides",
                                     "LFQ.intensity"))



# ---- Save data image ----
save.image("Data/RData/01.RData")
