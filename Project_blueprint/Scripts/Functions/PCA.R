#' Evaluates data cell-wise
#'
#' @param data_ data list
#' @param expr expression to be applied to each cell (must contain x as
#' variable, e.g. x > 2)
#' @param modify (optional) character vector of columns to modify
#' @param ignore (optional) character vector of columns to ignore
#' @param input.name if data_ is list: name of data to use
#' @param output.name if data_ is list: name of output data to save in list under
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#'
do_pca <- function(data_, 
                   input.name, 
                   output.name = "_pca") {
  
  # Check input
  data <- .unpack_data(data_, input.name)
  
  # Save attributes
  data_attributes <- attributes(data)
  
  ####
  
  
  # Compute PCA
  data.prcomp <- data %>%
    dplyr::select(dplyr::all_of(.data_columns(data, data_attributes))) %>%
    prcomp()
  
  data <- data %>%
    dplyr::select(-dplyr::all_of(.data_columns(data, data_attributes))) %>%
    cbind(data.prcomp[["x"]]) %>% 
    tibble::as_tibble()
  
  
  # Output name
  if (substr(output.name, 1, 1) == "_") {
    output.name <- paste0(data_attributes[["input.name"]], output.name)
  } 
  
  ####
  
  # Prepare return
  data_ <- .pack_data(data, data_, data_attributes, output.name, overwrite = F)
  
  # Return
  return(data_)
  
}