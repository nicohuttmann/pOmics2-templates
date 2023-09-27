

## Define folder structures
dir.create("Data")
dir.create("Data/RData")
dir.create("Scripts")
dir.create("Output")


# ---- Install packages ---- 
install.packages("tidyverse")
install.packages("devtools")
# to install the version we used during the workshop
devtools::install_github(repo = "nicohuttmann/pOmics2", 
                         ref = "02ed0624ce6a050ee5272c05d8fdec20701779d0")
# to install the most recent version
devtools::install_github("nicohuttmann/pOmics2")



install.packages("ggvenn")
install.packages("eulerr")

install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")

install.packages("openxlsx")
