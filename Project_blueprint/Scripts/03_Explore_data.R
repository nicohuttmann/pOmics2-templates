
#
# 03
# Overview analysis
# 



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Load data image ---- 
load("Data/RData/02.RData")



# ---- Define thresholds for identified proteins ----
# Create final lists of identified proteins in datasets
# and samples individually
Analysis[["Id_proteins"]] <- 
  get_data_frame(which = "Mean.peptides", 
                 variables = clean, 
                 dataset = "proteinGroups_c") %>% 
  do_expr(FUN = function(x) x > 0) %>% 
  add_observations_data("genotype") %>% 
  do_row_summary(function(x) sum(x) > 2, by = "genotype") %>% 
  do_transpose() %>% 
  do_fun(save_variables_data, 
         column = c("IRP1", 
                    "WT"), 
         name = c("ID_IRP1", 
                  "ID_WT"), 
         dataset = "proteinGroups_c")



# ---- Compare identified proteins ----
## with Venn diagrams
Analysis[["Compare_dataset_proteins"]] <- 
  get_variables_data(which = c("ID_WT", "ID_IRP1"), 
                     variables = ID_WT | ID_IRP1, 
                     dataset = "proteinGroups_c") %>% 
  plot_venn(c("ID_WT", "ID_IRP1"))

## with an Euler diagram
Analysis[["Compare_dataset_proteins"]] <- 
  get_variables_data(which = c("ID_WT", "ID_IRP1"), 
                     variables = ID_WT | ID_IRP1, 
                     dataset = "proteinGroups_c") %>% 
  plot_euler(c("ID_WT", "ID_IRP1"))


## Add gene and protein names to dataframe
Analysis[["Compare_dataset_proteins"]] <- 
  Analysis[["Compare_dataset_proteins"]] %>% 
  add_variables_data(which = "Gene.names", 
                     dataset = "proteinGroups_c", 
                     input.name = "variables_data") %>% 
  add_variables_data(which = "Protein.names", 
                     dataset = "proteinGroups_c")



# ---- Missing value imputation ---- 

## Count valid LFQ values per protein
## (this helps to know which proteins to keep for the analysis)
get_data_frame(which = "Mean.LFQ.intensity", 
               variables = clean, 
               dataset = "proteinGroups_c") %>% 
  do_row_summary(function(x) mean(x > 0)) %>% 
  do_transpose() %>% 
  do_fun(save_variables_data, 
         column = "all", 
         name = "nLFQ", 
         dataset = "proteinGroups_c")

## Impute missing values from normal distribution
## (there are many different ways to imputte values;
##  others are KNN, mean imputation, or Random Forest)
get_data_frame(which = "Mean.LFQ.intensity", 
               variables = nLFQ >= 0.5, 
               dataset = "proteinGroups_c") %>%  
  do_fun(impute_norm, 
         shift = 1.8, 
         width = 0.3) %>% 
  do_fun(save_data_frame, 
         name = "LFQ.imputed", 
         "proteinGroups_c")





# ---- PCA ---- 
Analysis[["PCA"]] <- 
  get_data_frame(which = "LFQ.imputed", 
                 variables = , 
                 dataset = "proteinGroups_c") %>% 
  do_expr(FUN = scale, eval.rowwise = F) %>% 
  do_pca() %>% 
  add_observations_data("genotype", dataset = "proteinGroups_c") %>% 
  plot_pca(x = "PC1", 
           y = "PC2", 
           color.by = "genotype", size = 2) %>% 
  plot_pca(x = "PC3", 
           y = "PC5", 
           color.by = "genotype", size = 2)


  # Analysis[["PCA"]] <- get_data_frame(which = "Mean.LFQ.intensity", 
  #                variables = clean & nLFQ == 1, 
  #                dataset = "proteinGroups_c") %>% 
  #   do_expr(FUN = scale, eval.rowwise = F) %>% 
  #   do_pca() %>% 
  #   add_observations_data("group", dataset = "proteinGroups_c") %>% 
  #   add_observations_data("genotype", dataset = "proteinGroups_c") %>% 
  #   add_observations_data("diet", dataset = "proteinGroups_c") %>% 
  #   plot_pca(x = "PC1", 
  #            y = "PC2", 
  #            color.by = "diet", 
  #            color = c("#4DAF4A", 
  #                      "#E41A1C",
  #                      "#377EB8"), 
  #            shape.by = "genotype", 
  #            shape = c(17, 16), size = 3) %>% 
  #   plot_pca(x = "PC3", 
  #            y = "PC5",
  #            color.by = "diet", 
  #            color = c("#4DAF4A", 
  #                      "#E41A1C",
  #                      "#377EB8"), 
  #            shape.by = "genotype", 
  #            shape = c(17, 16), size = 3)
    

## Save Principal component values for further analysis
Analysis[["PCA"]] %>% 
  do_fun(save_observations_data, 
         column = c("PC1", "PC2", "PC3", "PC4", "PC5"), 
         dataset = "proteinGroups_c", 
         input.name = "LFQ.imputed_expr_pca_obs")






# ---- t-Test/fold-change/Volcano plot ----
Analysis[["Volcano_WT-IRP1"]] <- 
  get_data_frame(which = "LFQ.imputed", 
                 variables = , 
                 observations = , 
                 dataset = "proteinGroups_c") %>% 
  add_observations_data("genotype") %>% 
  do_t.test(group.column = "genotype", 
            control.group = "WT", 
            fc.threshold = 1.5, 
            p.value.cutoff = 0.01, 
            paired = F, 
            var.equal = T, 
            pAdjustMethod = "BH") %>% 
  add_variables_data("Protein.names") %>% 
  add_variables_data("Gene.names") %>% 
  plot_volcano()




# ---- Heatmap ---- 
Analysis[["Heatmap"]] <- 
  get_data_frame(which = "LFQ.imputed", 
                 variables = , 
                 observations = , 
                 dataset = "proteinGroups_c") %>% 
  do_expr(FUN = scale, eval.rowwise = F) %>% 
  do_transpose() %>% 
  do_hclust_x() %>% 
  do_hclust_y(input.name = "LFQ.imputed_expr_t") %>% 
  plot_ComplexHeatmap(dend_x = "LFQ.imputed_expr_t_hclustx", 
                      dend_y = "LFQ.imputed_expr_t_hclusty", 
                      border = T, 
                      show_column_names = T, 
                      input.name = "LFQ.imputed_expr_t")



# ---- Save data image ---- 
save.image("Data/RData/03.RData")
