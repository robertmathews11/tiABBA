prepare_time <- function(input_expression){
  # Process individual cluster group
  gene_count <- nrow(input_expression)

  # Make data frame to return for plotting
  res <- data.frame(
    # Values for plotting
    expression = apply(input_expression, 2, mean),
    var = apply(input_expression, 2, var),
    time = colnames(input_expression) %>% str_split_i(pattern = "_", 2) %>% as.numeric()
  )
  res_extend <- res %>% filter(.data$time <= 8)
  res_extend <- res_extend %>% mutate(time = .data$time + 24)

  res <- rbind(res, res_extend)

  return(res)
}

time_plot <- function(input_time_data, ...){


  element_time_base(data = input_time_data, ...) +
    element_time_background(condition = unique(input_time_data$condition)) +
    element_time_ribbon(
      ymin = (input_time_data[["expression"]] - input_time_data[["var"]]),
      ymax = (input_time_data[["expression"]] + input_time_data[["var"]])
    ) +
    element_time_line() +
    element_time_breaks() +
    element_time_theme()
}
element_time_base <- function(data, x = data$time, y = data$expression, group = NULL, color = NULL,
                              linewidth = NULL, linetype = NULL) {
  ggplot(data, aes(x = {{ x }}, y = {{ y }}, group = {{ group }},
                   color = {{ color }}, fill = {{ color }},
                   linewidth = {{ linewidth }}, linetype = {{ linetype }}))
}
element_time_background <- function(condition = "none"){
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

element_time_ribbon <- function(ymin, ymax, set_alpha = 0.3){
  list(
    geom_ribbon(
      aes(ymin = pmax(ymin, -1), ymax = pmin(ymax, 1)),
      color = NA,
      alpha = set_alpha
    ),
    scale_fill_identity()
  )
}
element_time_line <- function(){
  list(
    geom_line(),
    scale_color_identity(),
    scale_linewidth_identity(),
    scale_linetype_identity()
  )
}

element_time_breaks <- function(){
  list(
    scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20,24,28,32), limits = c(0, 32)),
    scale_y_continuous(breaks = c(-1,0,1)),
    coord_cartesian(ylim = c(-1,1))
  )
}
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

make_violin_plot <- function(input_data, input_threshold, group_i) {

  good_above <- c(eigen = FALSE, gene = TRUE, sil = TRUE, var = FALSE)
  group_colors <- c(
    eigen = "#1B9E77",
    gene  = "#D95F02",
    sil   = "#7570B3",
    var   = "#E7298A"
  )

  dat <- input_data %>%
    filter(.data$group == group_i)

  threshold <- input_threshold %>%
    filter(.data$group == group_i) %>%
    pull(.data$threshold)

  dens <- density(dat$value, n = 512)
  dens_df <- data.frame(y = dens$x, width = dens$y / max(dens$y) * 0.45)

  # Interpolate width at the threshold so the color split has no gap
  width_at_threshold <- approx(dens_df$y, dens_df$width, xout = threshold)$y
  dens_df <- bind_rows(dens_df, data.frame(y = threshold, width = width_at_threshold)) %>%
    arrange(.data$y) %>%
    mutate(is_good = if (good_above[[group_i]]) .data$y >= threshold else .data$y <= threshold)

  good_df <- dens_df %>%
    filter(.data$is_good)
  bad_df  <- dens_df %>%
    filter(!.data$is_good)

  find_above <- nrow(bad_df) / (nrow(good_df) + nrow(bad_df)) * 100

  find_above <- ifelse(find_above > 1, round(find_above, 0),
                       99)

  set_title <- str_c(str_to_title(group_i), "\n\n",
                     find_above, "% Rank", "\n\n",
                     round(threshold, 3))

  ggplot() +
    geom_ribbon(data = bad_df, aes(y = bad_df$y, xmin = -bad_df$width, xmax = bad_df$width),
                orientation = "y", fill = "#CCCCCC") +
    geom_ribbon(data = good_df, aes(y = good_df$y, xmin = -good_df$width, xmax = good_df$width),
                orientation = "y", fill = group_colors[[group_i]]) +
    geom_hline(yintercept = threshold, linetype = 3, linewidth = 1) +
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    labs(y = set_title) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = rel(1.1), angle = 0, vjust = 0.5),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = rel(1.1)),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}

find_iqr <- function(input_vector){
  iqr <- IQR(input_vector)
  q1 <- quantile(input_vector, 0.25)
  q3 <- quantile(input_vector, 0.75)

  low_lim <- (q1 - (1.5 * iqr))
  up_lim <- (q3 + (1.5 * iqr))

  return(c(low_lim, up_lim))
}

