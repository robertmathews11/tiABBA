#' Violin plot creation
#'
#' This function creates a violin plot based on relative rankings of the parameter scan
#'
#' @param input_data A data.frame object contained indexed results from [tiABBA::scan_index]
#' @param input_threshold The value to be compared against other results, such as percent genes retained
#' @param group_i The metric being compared
#'
#' @returns A geom_violin plot
#'
#' @export

make_violin_plot <- function(

  input_data,
  input_threshold,
  group_i

) {

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

  res <- ggplot() +
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

  suppressWarnings(print(res))
}
