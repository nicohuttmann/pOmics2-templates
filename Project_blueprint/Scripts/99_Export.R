
#
# 04
# Comparative analysis
# 



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)
library(openxlsx)


# ---- Load data image ---- 
openxlsx::write.xlsx(x = Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var"]], 
                     file = "Output/Volcano_plot_data.xlsx")

openxlsx::write.xlsx(x = list(Volcano = Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var"]] %>% 
                                dplyr::arrange(p.adjust), 
                              Upregulated = Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var"]] %>% 
                                dplyr::filter(regulated == "up") %>% 
                                dplyr::arrange(desc(log2.fc)), 
                              Downregulated = Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var"]] %>% 
                                dplyr::filter(regulated == "down") %>% 
                                dplyr::arrange(log2.fc), 
                              ORA_up = data.frame(Analysis[["Volcano_WT-IRP1"]][["LFQ.imputed_obs_t.test_var_var_with"]])), 
                     file = "Output/Comparative_abundance_analysis.xlsx")



# ---- Export figures ---- 

## You can simply navigate to the 'Plots' panel to your right
## and manually click 'Export' or as described below

## PNG
png(filename = "Output/PCA_1_2.png")
# Plot your graphic
Analysis[["PCA"]][["LFQ.imputed_expr_pca_obs_plotpca"]]
# Close the graphics device
dev.off()


## With pOmics2 package function
export_pdf(p = Analysis[["Volcano_WT-IRP1"]][["plot_volcano"]], 
           file = "Output/Volcano_plot.pdf")



# ---- Save data image ---- 
#save.image("Data/RData/99.RData")
