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
do_t.test <-function(data_, 
                     group.column,
                     control.group,
                     paired = F,
                     var.equal = T,
                     p.value.cutoff = 0.05,
                     fc.threshold = 0, 
                     pAdjustMethod = "BH",
                     input.name, 
                     output.name = "_t.test") {
  
  # Check input
  data <- .unpack_data(data_)
  
  # Save attributes
  data_attributes <- attributes(data)
  
  ####
  
  # Prepare results data frame
  data_t.test <- dplyr::tibble(variables = .data_columns(data, data_attributes),
                               log2.fc = NA_real_,
                               sig.log2.fc = NA,
                               p.value = NA_real_,
                               sig.p.value = NA,
                               significant = NA,
                               regulated = "not")
  
  # Check group information
  if (!hasArg(group.column) || 
      (length(group.column) != 1) || 
      !is.character(group.column)) {
    stop("The <group.column> must be a character indicating which column", 
         " contains information of group membership for each observation.", 
         call. = FALSE)
  }
  
  # Group information
  group_info <- data[[group.column]]
  
  if (length(unique(group_info)) != 2) {
    stop("The <group.column> column must contain only two different groups.", 
         call. = FALSE)
  }
  
  
  if (!hasArg(control.group)) {
    stop("Please specify the <control.group> of your comparison.", 
         call. = FALSE)
  } else if (!control.group %in% group_info) {
    stop("The <group.column> column must contain the <control.group>.", 
         call. = FALSE)
  }
  
  group_info <- factor(group_info, 
                       levels = c(control.group, 
                                  setdiff(unique(group_info), control.group)))
  
  
  list_t.test <- list()
  
  
  for (i in data_t.test[["variables"]]) {
    
    # t-test model
    list_t.test[[i]] <- t.test(log2(data[[i]][group_info == levels(group_info)[2]]),
                               log2(data[[i]][group_info == levels(group_info)[1]]),
                               var.equal = var.equal,
                               paired = paired)
    
    # log2 fold-change
    data_t.test[match(i, data_t.test[["variables"]]), "log2.fc"] <-
      log2(mean(data[[i]][group_info == levels(group_info)[2]]) /
             mean(data[[i]][group_info == levels(group_info)[1]]))
    
    # p-value
    data_t.test[match(i, data_t.test[["variables"]]), "p.value"] <-
      list_t.test[[i]]$p.value
    
  }
  
  
  ## Determine significant effects
  
  # Numeric threshold
  if (is.numeric(fc.threshold)) {
    
    # Absolute value given
    if (length(fc.threshold) == 1) {
      
      data_t.test[, "sig.log2.fc"] <- abs(data_t.test[, "log2.fc"]) > fc.threshold
      
      # Values for positive and negative threshold given
    } else {
      
      data_t.test[, "sig.log2.fc"] <-
        ifelse(data_t.test[, "log2.fc"] > 0,
               data_t.test[, "log2.fc"] > fc.threshold[2],
               data_t.test[, "log2.fc"] < fc.threshold[1])
      
    }
    
    
    # Determine significant effect size by confidence interval
  } else if (is.character(fc.threshold) && grep(fc.threshold, "confidence|interval")) {
    
    for (i in data_t.test[["variables"]]) {
      data_t.test[match(i, data_t.test[["variables"]]), "sig.log2.fc"] <-
        sum(list_t.test[[i]]$conf.int > 0) != 1
    }
    
  } else {
    stop('The given <fc.threshold> is not supported. Either provide a numeric threshold (e.g. 1), a vector containing a negative and a positive" threshold (e.g. c(-1, 2)) or indicate "interval".', 
         call. = FALSE)
  }
  
  # p-value
  if (is.character(p.value.cutoff)) stop("Please provide a numeric p-value cutoff.")
  
  # Adjust p-value
  data_t.test <- data_t.test %>%
    # Adjust p-value
    dplyr::mutate(p.adjust = p.adjust(p.value, method = pAdjustMethod),
                  .after = p.value) %>%
    # Significant p-values
    dplyr::mutate(sig.p.value = p.adjust < p.value.cutoff,
                  .after = p.adjust) %>%
    # Final significance
    dplyr::mutate(significant = sig.log2.fc & sig.p.value) %>%
    #
    dplyr::rowwise() %>%
    dplyr::mutate(regulated = if (significant) ifelse(log2.fc > 0, "up", "down")
                  else "not") %>%
    dplyr::ungroup()
  
  
  
  # Output name
  if (substr(output.name, 1, 1) == "_") {
    output.name <- paste0(data_attributes[["input.name"]], output.name)
  } 
  
  ####
  
  # Prepare return
  data_ <- .pack_data(data_t.test, data_, data_attributes, output.name, overwrite = F)
  
  # Return
  return(data_)
  
}


#' Plots a volcano plot with ggplot2
#'
#' @param data_ data_
#' @param color column for color of points
#' @param p.value.cutoff p-value limit for coloring
#' @param pos.log2fc.cutoff positive log2 fold-change limit for coloring
#' @param neg.log2fc.cutoff negative log2 fold-change limit for coloring
#' @param highlight.variables variables to highlight by point.size
#' @param highlight.color color to use to highlight proteins
#' @param x.axis.title title of x-axis
#' @param y.axis.title title of y-axis
#' @param text.size size of text in points (5-8)
#' @param text.color color of text
#' @param point.size point size (0.5-2)
#' @param point.alpha transparency (0-1)
#' @param highlight.point.size size of point of highlighted variables
#' @param highlight.point.alpha transparency of points to be highlighted
#' @param x.axis.breaks break size between ticks of x-axis
#' @param y.axis.breaks break size between ticks of y-axis
#' @param axis.line.size width of axes lines
#' @param axis.color color of axes lines
#' @param axis.ticks.size width of axis ticks
#' @param axis.title.size size of axis title
#' @param axis.text.size size of axis labels
#' @param aspect.ratio y/x ratio
#' @param use.plotly make interactive plots with ggplotly() (default: T if html
#' document)
#' @param view view plot
#' @param input name of input data
#' @param output name of output data
#'
#' @return
#' @export
#'
#'
#'
plot_volcano <- function(data_,
                         color = "regulated",
                         p.value.cutoff = 0.05,
                         pos.log2fc.cutoff = 0,
                         neg.log2fc.cutoff = 0,
                         highlight.variables = NULL,
                         highlight.color = NULL,
                         x.axis.title = "log2 fold-change",
                         y.axis.title = "-log10(p-value)",
                         text.size = 6,
                         text.color = "black",
                         point.size = 2,
                         point.alpha = 0.8,
                         highlight.point.size = 3,
                         highlight.point.alpha = 0.8,
                         x.axis.breaks = 1,
                         y.axis.breaks = 1,
                         axis.line.size = 0.5,
                         axis.color = "black",
                         axis.ticks.size = 0.3,
                         axis.title.size = 8,
                         axis.text.size = 6,
                         aspect.ratio = 0.8,
                         use.plotly = F,
                         input.name, 
                         output.name = "plot_volcano") {
  
  # Check input
  data <- .unpack_data(data_, input.name)
  
  # Save attributes
  data_attributes <- attributes(data)
  
  ####
  
  
  # Add color column
  data <- data %>%
    dplyr::rowwise() %>%
    dplyr::rename(color = !!color) %>%
    {if (!is.null(highlight.color))
      dplyr::mutate(., color = ifelse(variables %in% highlight.variables,
                                      "highlight",
                                      color))
      else .} %>%
    dplyr::mutate(point.size =
                    ifelse(variables %in% highlight.variables,
                           highlight.point.size,
                           point.size)) %>%
    dplyr::mutate(alpha =
                    ifelse(variables %in% highlight.variables,
                           highlight.point.alpha,
                           point.alpha)) %>%
    dplyr::arrange(point.size, -p.value)
  
  
  
  # Add text column for interactive plots
  if (use.plotly) {
    data$text <- paste(p2g(data$variables),
                       p2n(data$variables),
                       data$variables, sep = "\n")
  } else {
    #data$text <- p2g(data$variables)
    data$text <- "" #data$variables
  }
  
  
  # Line width factor
  lwf <- 1 / (ggplot2::.pt * 72.27 / 96)
  
  
  p <- ggplot(data = data,
              aes(x = log2.fc,
                  y = -log10(p.value),
                  col = color,
                  size = point.size,
                  text = text,
                  alpha = alpha)) +
    geom_point(shape = 16, stroke = 0) +
    scale_size_continuous(range = range(data$point.size)) +
    scale_alpha_identity() +
    theme(aspect.ratio = aspect.ratio,
          axis.line = element_line(size = 0.5 * lwf, color = axis.color),
          panel.background = element_blank(),
          axis.text = element_text(size = axis.text.size, color = text.color),
          axis.ticks = element_line(size = axis.ticks.size * lwf,
                                    color = axis.color),
          axis.title = element_text(size = axis.title.size, color = text.color),
          legend.position = 0) +
    scale_color_manual(values = c("highlight" = highlight.color,
                                  "up" = "red",
                                  "down" = "blue",
                                  "not" = "grey")) +
    scale_x_continuous(
      limits = .axis_limit_breaks(plot.limits = range(data$log2.fc),
                                  break.size = x.axis.breaks)$limits,
      breaks = .axis_limit_breaks(plot.limits = range(data$log2.fc),
                                  break.size = x.axis.breaks)$breaks,
      expand = c(0, 0)) +
    scale_y_continuous(
      limits = .axis_limit_breaks(plot.limits = range(-log10(data$p.value)),
                                  break.size = y.axis.breaks)$limits,
      breaks = .axis_limit_breaks(plot.limits = range(-log10(data$p.value)),
                                  break.size = y.axis.breaks)$breaks,
      expand = c(0, 0)) +
    xlab(x.axis.title) +
    ylab(y.axis.title)
  
  
  # Generate interactive plot with ggplotly()
  if (use.plotly) {
    p <- plotly::ggplotly(p, tooltip = "text")
  }
  
  # Print to Plots panel
  plot(p)
  
  # Output name
  if (substr(output.name, 1, 1) == "_") {
    output.name <- paste0(data_attributes[["input.name"]], output.name)
  }
  
  ####
  
  # Prepare return
  data_ <- .pack_data(p, data_, data_attributes, output.name, 
                      overwrite = T, 
                      last.entry = data_attributes[["input.name"]])
  
  # Return
  return(invisible(data_))
  
}
