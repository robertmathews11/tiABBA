#' Plot time data
#'
#' @param input_time_data A data.frame object prepared with [tiABBA::prepare_time]
#' @param ... Additional parameters passed on to ggplot elements
#'
#' @returns One png line plot for group-mean expression
#'
#' @export

time_plot <- function(

  input_time_data,

  ...

) {


  res <- element_time_base(data = input_time_data, ...) +
    element_time_background(condition = unique(input_time_data$condition)) +
    element_time_ribbon(
      ymin = (input_time_data[["expression"]] - input_time_data[["var"]]),
      ymax = (input_time_data[["expression"]] + input_time_data[["var"]])
    ) +
    element_time_line() +
    element_time_breaks() +
    element_time_theme()

  suppressWarnings(print(res))
}


#' Element plot base
#'
#' @param data A data.frame input to build the basic ggplot from
#' @param x The x aes in ggplot
#' @param y The y aes in ggplot
#' @param group Designates the column which identifies group to connect the geom_line
#' @param color Designates the column which identifies color in scale_color_identity
#' @param linewidth Designates the column which identifies linewidth in scale_linewidth_identity
#' @param linetype Designates the column which identifies linetype in scale_linetype_identity
#'
#' @returns A ggplot frame
#'
#' @export

element_time_base <- function(

  data,
  x = data$time,
  y = data$expression,
  group = NULL,
  color = NULL,
  linewidth = NULL,
  linetype = NULL

) {
  ggplot(data, aes(x = {{ x }}, y = {{ y }}, group = {{ group }},
                   color = {{ color }}, fill = {{ color }},
                   linewidth = {{ linewidth }}, linetype = {{ linetype }}))
}


#' Element plot background
#'
#' @param condition Optional parameter used to tune condition-specific contextual backgrounds
#'
#' @returns Adds a background to the [tiABBA::element_time_base] frame
#'
#' @export

element_time_background <- function(

  condition = "none"

) {
  list(
    annotate("rect", xmin = 12, xmax = 24, ymin = -1, ymax = 1, fill = "gray50", alpha = 0.5),
    annotate("rect", xmin = 0, xmax = 12, ymin = -1, ymax = 1, fill = "gray90", alpha = 0.5),
    annotate("rect", xmin = 24, xmax = 32, ymin = -1, ymax = 1, fill = "gray90", alpha = 0.5),
    {if("tan" %in% condition) annotate("rect", xmin = 13, xmax = 15, ymin = -1, ymax = 1, fill = "red", alpha = 0.3)},
    {if("tan" %in% condition) annotate("rect", xmin = 21, xmax = 23, ymin = -1, ymax = 1, fill = "red", alpha = 0.3)},
    annotate("segment", x = 12, y = -1, yend = 1, linetype = 3, linewidth = 1.5, color = "black"),
    annotate("segment", x = 24, y = -1, yend = 1, linetype = 3, linewidth = 1.5, color = "black"),
    annotate("segment", x = 0, y = 0, xend = 32, linetype = 2, linewidth = 1.5, color = "black")
  )
}


#' Element geom_ribbon
#'
#' @param ymin Lower boundary of the geom_ribbon
#' @param ymax Upper boundary of the geom_ribbon
#' @param set_alpha transparency of the geom_ribbon
#'
#' @returns Adds a geom_ribbon to the [tiABBA::element_time_base] frame
#'
#' @export

element_time_ribbon <- function(

  ymin,
  ymax,
  set_alpha = 0.3

) {
  list(
    geom_ribbon(
      aes(ymin = pmax(ymin, -1), ymax = pmin(ymax, 1)),
      color = NA,
      alpha = set_alpha
    ),
    scale_fill_identity()
  )
}


#' Element geom_line
#'
#' @returns Adds a geom_line to the [tiABBA::element_time_base] frame
#'
#' @export

element_time_line <- function() {
  list(
    geom_line(),
    scale_color_identity(),
    scale_linewidth_identity(),
    scale_linetype_identity()
  )
}


#' Element axis scales
#'
#' @returns Adjusts the axis scales and breaks in the [tiABBA::element_time_base] frame
#'
#' @export

element_time_breaks <- function() {
  list(
    scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20,24,28,32), limits = c(0, 32)),
    scale_y_continuous(breaks = c(-1,0,1)),
    coord_cartesian(ylim = c(-1,1))
  )
}


#' Element plot theme
#'
#' @param x_lab Label to be used in the x axis
#' @param y_lab Label to be used in the y axis
#'
#' @returns Adjusts the theme and axis labels in the [tiABBA::element_time_base] frame
#'
#' @export

element_time_theme <- function(x_lab = NULL, y_lab = NULL){
  list(
    theme_minimal(),
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10)),
    labs(x = x_lab,
         y = y_lab)
  )
}
