
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
import_files(files = c("Data/raw/protein_PTM.tsv",
                       "Data/raw/protein_SP2.tsv"))



# ---- Tidy up data sets ----

## Regular preparation
Datasets[["SP2"]] <-
  import2new_dataset(raw.data = Info$Imports[["protein_SP2"]],
                     variable.identifiers = Protein.ID,
                     variables.data = c("Protein",
                                        "Protein.ID",
                                        "Entry.Name",
                                        "Gene",
                                        "Length",
                                        "Organism",
                                        "Protein.Description",
                                        "Protein.Existence",
                                        "Coverage",
                                        "Protein.Probability",
                                        "Top.Peptide.Probability",
                                        "Total.Peptides",
                                        "Unique.Peptides",
                                        "Razor.Peptides",
                                        "Total.Spectral.Count",
                                        "Unique.Spectral.Count",
                                        "Razor.Spectral.Count",
                                        "Total.Intensity",
                                        "Unique.Intensity",
                                        "Razor.Intensity",
                                        "Razor.Assigned.Modifications",
                                        "Razor.Observed.Modifications",
                                        "Indistinguishable.Proteins"),
                     observations = c("brain_total",
                                      "eyes_total",
                                      "heart_total",
                                      "liver_total",
                                      "muscle_total",
                                      "brain_soluble",
                                      "eyes_soluble",
                                      "heart_soluble",
                                      "liver_soluble",
                                      "muscle_soluble"),
                     data.frames = c("TMT"))

## PTM preparation
Datasets[["PTM"]] <-
  import2new_dataset(raw.data = Info$Imports[["protein_PTM"]],
                     variable.identifiers = Protein.ID,
                     variables.data = c("Protein",
                                        "Protein.ID",
                                        "Entry.Name",
                                        "Gene",
                                        "Length",
                                        "Organism",
                                        "Protein.Description",
                                        "Protein.Existence",
                                        "Coverage",
                                        "Protein.Probability",
                                        "Top.Peptide.Probability",
                                        "Total.Peptides",
                                        "Unique.Peptides",
                                        "Razor.Peptides",
                                        "Total.Spectral.Count",
                                        "Unique.Spectral.Count",
                                        "Razor.Spectral.Count",
                                        "Total.Intensity",
                                        "Unique.Intensity",
                                        "Razor.Intensity",
                                        "Razor.Assigned.Modifications",
                                        "Razor.Observed.Modifications",
                                        "Indistinguishable.Proteins"),
                     observations = c("brain_total",
                                      "eyes_total",
                                      "heart_total",
                                      "liver_total",
                                      "muscle_total",
                                      "brain_soluble",
                                      "eyes_soluble",
                                      "heart_soluble",
                                      "liver_soluble",
                                      "muscle_soluble"),
                     data.frames = c("TMT"))



# ---- Save data image ----
save.image("Data/RData/01.RData")