
#
# 03
# Overview analysis
#



# ---- Load libraries ----
library(pOmics2)
library(tidyverse)



# ---- Load data image ----
load("Data/RData/02.RData")



# ---- Compare identified proteins ---- 

## IDs
get_variables_data(which = "min2",  
                   dataset = "SP2") %>% 
  get_variables_data(which = "min2", 
                     dataset = "PTM") %>% 
  do_join(join.type = "full", by = "variables", 
          input.names = c(1, 2)) %>% 
  do_expr(FUN = function(x) ifelse(is.na(x), F, x)) %>% 
  plot_euler(columns = c(SP2 = "min2_1", PTM = "min2_2"))


## By total intensity in datasets
get_variables_data(which = "Unique.Intensity",  
                   dataset = "SP2") %>% 
  do_expr(FUN = log10) %>% 
  get_variables_data(which = "Unique.Intensity", 
                     dataset = "PTM") %>% 
  do_expr(FUN = log10) %>% 
  do_join(join.type = "full", by = "variables", 
          input.names = c(2, 4)) %>% 
  plot_xy(x = "Unique.Intensity_2", 
          y = "Unique.Intensity_4")



# ----  Clean TMT data ---- 

## SP2
get_data_frame(which = "TMT", dataset = "SP2") %>% 
  do_transpose() %>% 
  do_column_summary(function(x) mean(x > 0), name = "TMT_f") %>% 
  do_fun(save_variables_data, 
         column = "TMT_f", 
         dataset = "SP2")

# View
get_variables_data(which = "TMT_f", 
                   dataset = "SP2") %>% 
  do_with(table(TMT_f))

# Threshold protein quantifications and save to dataset
get_variables_data(which = "TMT_f",  
                   dataset = "SP2") %>% 
  do_expr(FUN = function(x) x == 1) %>% 
  do_with(table(TMT_f)) %>% 
  do_fun(save_variables_data, 
         column = "TMT_f", 
         name = "TMT0", 
         dataset = "SP2", 
         fill = NA, 
         input.name = 2) 


## PTM
get_data_frame(which = "TMT", dataset = "PTM") %>% 
  do_transpose() %>% 
  do_column_summary(function(x) mean(x > 0), name = "TMT_f") %>% 
  do_fun(save_variables_data, 
         column = "TMT_f", 
         dataset = "PTM")

get_variables_data(which = "TMT_f", 
                   dataset = "PTM") %>% 
  do_with(table(TMT_f))

# Threshold protein quantifications and save to dataset
get_variables_data(which = "TMT_f",  
                   dataset = "PTM") %>% 
  do_expr(FUN = function(x) x == 1) %>% 
  do_with(table(TMT_f)) %>% 
  do_fun(save_variables_data, 
         column = "TMT_f", 
         name = "TMT0", 
         dataset = "PTM", 
         fill = NA, 
         input.name = 2)



# ---- Heatmaps ---- 
get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  do_expr(FUN = log10) %>% 
  do_transpose() %>% 
  plot_ComplexHeatmap()

get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "PTM") %>% 
  do_expr(FUN = log10) %>% 
  do_transpose() %>% 
  plot_ComplexHeatmap()



# ---- Boxplot ---- 


data_ <- get_data_frame(which = "TMT", 
                        variables = TMT0, 
                        dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  add_observations_data("tissue") %>% 
  plot_gg(aes(x = tissue, y = log10(value), fill = fraction), 
          geom_boxplot()) %>% 
  get_data_frame(which = "TMT", 
                 variables = TMT0, 
                 dataset = "PTM") %>% 
  add_observations_data("fraction") %>% 
  add_observations_data("tissue") %>% 
  plot_gg(aes(x = tissue, y = log10(value), fill = fraction), 
          geom_boxplot()) %>% 
  plot_patchwork(input.names = c(4, 8))


get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun_grouped("fraction", 
                 norm_vsn, 
                 input.name = 2, 
                 output.name = "vsn separate") %>% 
  add_observations_data("tissue") %>% 
  do_expr(FUN = log10) %>% 
  plot_gg(aes(x = fraction, y = value, fill = fraction), 
          geom_boxplot(), 
          facet_grid(cols = vars(tissue)))





get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun_grouped("fraction", 
                 norm_vsn, 
                 input.name = 2, 
                 output.name = "vsn separate") %>% 
  add_observations_data("tissue") %>% 
  do_expr(FUN = log10) %>% 
  plot_gg(aes(x = tissue, y = value, fill = fraction), 
          geom_boxplot())






a <- data %>% 
  column_to_rownames(var = "observations") %>% 
  as.matrix() %>% 
  log2() %>% 
  t() %>% 
  vsn::meanSdPlot()

b <- vsn::meanSdPlot(eset, plot = F)

b$gg



get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  #add_observations_data("fraction") %>% 
  do_fun(norm_vsn) %>% 
  do_expr(FUN = log10) %>% 
  last() %>% 
  pivot_longer(cols = -1) %>% 
  ggplot(aes(x = observations, y = value)) + 
  geom_boxplot()

get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  #add_observations_data("fraction") %>% 
  do_fun(norm_limma_MedianValues) %>% 
  do_expr(FUN = log10) %>% 
  last() %>% 
  pivot_longer(cols = c(-1, -2)) %>% 
  ggplot(aes(x = observations, y = value)) + 
  geom_boxplot()


data <- get_data_frame(which = "TMT", 
                       variables = TMT0, 
                       dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun(norm_vsn_grouped, group.column = "fraction") %>% 
  do_expr(FUN = log10) %>% 
  do_transpose() %>% 
  plot_ComplexHeatmap()

data <- get_data_frame(which = "TMT", 
                       variables = TMT0, 
                       dataset = "PTM") %>% 
  add_observations_data("fraction") %>% 
  do_fun(norm_vsn_grouped, group.column = "fraction") %>% 
  do_expr(FUN = log10) %>% 
  do_transpose() %>% 
  plot_ComplexHeatmap()





get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun(norm_vsn, output.name = "vsn together") %>% 
  do_fun_grouped("fraction", 
                 norm_vsn, 
                 input.name = 2, 
                 output.name = "vsn separate") %>% 
  plot_vsn_meanSdPlot_m(log2 = T, 
                        input.names = c(raw = 2, tog = 3, sep = 4)) %>% 
  do_select(-"fraction", input.name = 4)




get_data_frame(which = "TMT", 
               variables = TMT0, 
               dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_select(-fraction, select.data.only = F) %>% 
  do_fun(save_data_frame(name = "TMT_vsn", 
                         dataset = "SP2"))









# ---- Boxplot ---- 
data <- get_data_frame(which = "TMT", 
                       variables = TMT0, 
                       output.type = "tibble", 
                       dataset = "SP2") %>% 
  do_expr(log10(x)) %>% 
  do_mutate("observations", factor(observations, ))

data %>% 
  pivot_longer(cols = -1) %>% 
  ggplot(aes(x = observations, y = value)) + 
  geom_boxplot()


data <- get_data_frame(which = "TMT", 
                       variables = TMT0, 
                       dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun_grouped("fraction", 
                 norm_vsn, 
                 input.name = 2, 
                 output.name = "vsn separate") %>% 
  add_observations_data("tissue") %>% 
  do_expr(FUN = log10) %>% 
  last()

data %>% 
  pivot_longer(cols = -c(1:3)) %>% 
  ggplot(aes(x = tissue, y = value, fill = fraction)) + 
  geom_boxplot() +
  theme_bw()






data_ <- get_data_frame(which = "TMT", 
                        variables = TMT0, 
                        dataset = "SP2") %>% 
  add_observations_data("fraction") %>% 
  do_fun_grouped("fraction", 
                 norm_limma_normalizeBetweenArrays, 
                 input.name = 2, 
                 output.name = "vsn separate") %>% 
  add_observations_data("tissue") %>% 
  do_expr(FUN = log10) %>% 
  plot_gg(aes(x = tissue, y = value, fill = factor(fraction)), 
          geom_boxplot())




a <- get_data_frame(which = "TMT", 
                    variables = TMT0, 
                    output.type = "matrix", 
                    dataset = "SP2") %>% 
  t()

b <- get_data_frame(which = "TMT", 
                    variables = TMT0, 
                    output.type = "matrix", 
                    dataset = "SP2") %>% 
  t() %>% 
  limma::normalizeBetweenArrays(method = "scale")




# ---- Save data image ----
save.image("Data/RData/03.RData")

