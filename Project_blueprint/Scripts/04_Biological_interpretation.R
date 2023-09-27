
#
# 04
# Comparative analysis
# 



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)
library(org.Mm.eg.db)



# ---- Load data image ---- 
load("Data/RData/03.RData")



# ---- Over-representation analysis (of upregulated proteins) ---- 
Analysis[["Volcano_WT-IRP1"]] <- 
  Analysis[["Volcano_WT-IRP1"]] %>% 
  do_with(clusterProfiler::enrichGO(gene = variables[regulated == "up"], 
                                    OrgDb = "org.Mm.eg.db", 
                                    keyType = "UNIPROT", 
                                    universe = variables), 
          input.name = "LFQ.imputed_obs_t.test_var_var")

## View the results
Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var_with"]] %>% 
  data.frame() %>% 
  View()

## Represent data as dotplot
Analysis[["Volcano_WT-IRP1"]] %>% 
  do_fun(clusterProfiler::dotplot)



# ---- Gene set enrichment analysis ----
Analysis[["Volcano_WT-IRP1"]] <- 
  Analysis[["Volcano_WT-IRP1"]] %>% 
  do_with(clusterProfiler::gseGO(geneList = sort(rank(na.omit(setNames(-log2.fc, variables))), 
                                                 decreasing = T), 
                                 OrgDb = "org.Mm.eg.db", 
                                 keyType = "UNIPROT", 
                                 scoreType = "pos"), 
          input.name = "LFQ.imputed_obs_t.test_var_var", 
          output.name = "GSEA_log2fc_up")

## Represent data as dotplot
Analysis[["Volcano_WT-IRP1"]] %>% 
  do_fun(clusterProfiler::dotplot)


# ---- Gene set enrichment analysis from PCA values ----
## Extract values (Principal components) from PCA results
Analysis[["PCA_cor"]] <- get_data_frame(which = "LFQ.imputed", 
                                        variables = ,
                                        dataset = "proteinGroups_c") %>% 
  add_observations_data(c("PC1")) %>% 
  add_observations_data(c("PC2")) %>% 
  do_cor(x = "PC1", output.name = "cor_PC1") %>% 
  add_variables_data(which = "Protein.names")

## Do the Gene set enrichment analysis
## (This function can be adjusted to other values
##  you may produce from statistical tests like p-values, 
##  or the general abundance of proteins)
Analysis[["PCA_cor"]] <- Analysis[["PCA_cor"]] %>% 
  do_with(clusterProfiler::gseGO(geneList = sort(rank(na.omit(setNames(cor, variables))), 
                                                 decreasing = T), 
                                 OrgDb = "org.Mm.eg.db", 
                                 keyType = "UNIPROT", 
                                 scoreType = "pos"), 
          input.name = "cor_PC1_var", 
          output.name = "PC1_pos")



# ---- Retrieve informations about proteins form databases ---- 
## Load library
library(AnnotationDbi)
## Load organism database
library(org.Mm.eg.db)
# You can find other databases here:
# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb
# Human database:
# #library(org.Hs.eg.db)


## Get information column names
columns(org.Mm.eg.db)

## Get possible keytypes to query the database
keytypes(org.Mm.eg.db)

## Get possible keys for an ID (example: UNIPROT)
keys(x = org.Mm.eg.db, 
     keytype = "UNIPROT")

## Query the database
Analysis[["Protein_information"]] <- 
  AnnotationDbi::select(x = org.Mm.eg.db, 
                        keys = get_variables(dataset = "proteinGroups_c")[1:50], 
                        columns = c("SYMBOL", 
                                    "GENENAME"), 
                        keytype = "UNIPROT")

View(Analysis[["Protein_information"]])

## Gene ontology annotations for proteins
## (this table will have redundant proteins entries, 
##  as each proteins may have multiple annotations)
Analysis[["Gene_ontology_information"]] <- 
  AnnotationDbi::select(x = org.Mm.eg.db, 
                        keys = get_variables(dataset = "proteinGroups_c")[1:50], 
                        columns = c("GOALL"), 
                        keytype = "UNIPROT")

View(Analysis[["Gene_ontology_information"]])

## Other databases work the same way. Examples are
## UniProt: https://bioconductor.org/packages/release/bioc/html/UniProt.ws.html
## Gene Ontology: https://bioconductor.org/packages/release/data/annotation/html/GO.db.html
  


# ---- Save data image ---- 
save.image("Data/RData/04.RData")
