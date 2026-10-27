#' Prepare directories
#'
#' This helper function prepares directories for saving data
#'
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#'
#' @return This function creates file directories and has no return
#'
#' @export

prepare_directory <- function(

  group_file

  ) {

  if(!dir.exists("app")){dir.create("app")}
  if(!dir.exists("app/data")){dir.create("app/data")}
  if(!dir.exists("app/images")){dir.create("app/images")}

  outdir <- file.path("app/data", group_file)
  filedir <- file.path("app/images", group_file)
  if(!dir.exists(outdir)){dir.create(outdir)}
  if(!dir.exists(filedir)){dir.create(filedir)}

  if(!dir.exists(file.path(filedir,"heatmap"))){dir.create(file.path(filedir,"heatmap"))}
  if(!dir.exists(file.path(filedir,"sil_k"))){dir.create(file.path(filedir,"sil_k"))}
  if(!dir.exists(file.path(filedir,"gap_0p5"))){dir.create(file.path(filedir,"gap_0p5"))}
  if(!dir.exists(file.path(filedir,"gap_0p6"))){dir.create(file.path(filedir,"gap_0p6"))}
  if(!dir.exists(file.path(filedir,"gap_0p7"))){dir.create(file.path(filedir,"gap_0p7"))}
  if(!dir.exists(file.path(filedir,"gap_0p8"))){dir.create(file.path(filedir,"gap_0p8"))}
  if(!dir.exists(file.path(filedir,"gap_0p9"))){dir.create(file.path(filedir,"gap_0p9"))}

  return("Directory setup complete")
}


#' Prepare expression data
#'
#' This helper function averages sample data by time point for use in expression plots
#'
#' @param expr A data frame or matrix object, with rows of gene expression readings and columns of sample-time names.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#'
#' @return A data frame matching the input row length, and columns averaged by time point
#'
#' @export

prepare_expr <- function(

  expr,
  group_file

  ) {

  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  res <- normalize_rows(expr)

  # Save data
  saveRDS(res, file.path("app/data", group_file, "expression"))

  return(paste("Finished Time Normalization, saved to", file.path("app/data",group_file,"expression")))
}


#' Min-Max normalize data
#'
#' This helper function applied min-max normalization to the input data
#'
#' @param expr A data frame or matrix object to be normalized
#'
#' @return A matrix with rows min-max normalized between -1 and 1
#'
#' @export

normalize_rows <- function(

  expr

  ) {
  # Apply normalization to each row
  normalized_expr <- t(apply(expr, 1, function(row) {
    min_val <- min(row)
    max_val <- max(row)

    # Normalize the row values between -1 and 1
    normalized_row <- (row - min_val) / (max_val - min_val) * 2 - 1
    return(normalized_row)
  }))

  # Convert the result back to a data frame
  normalized_expr <- as.matrix(normalized_expr)

  return(normalized_expr)
}


#' Make parameter file ID
#'
#' @param cor_cut Numeric value of correlation parameter
#' @param sil_cut Numeric value of silhouette parameter
#'
#' @returns Transformed, formatted parameter ID
#'
#' @export

make_param_id <- function(

  cor_cut,
  sil_cut

  ) {
  fmt <- function(x) {
    x <- formatC(as.numeric(x), format = "f", digits = 2)
    x <- gsub("-", "n", x)
    gsub("\\.", "p", x)
  }
  paste0("cor", fmt(cor_cut), "_sil", fmt(sil_cut))
}


#' Downsample large plot data
#'
#' @param sub_mat A subsetting matrix
#' @param max_dim Maximum dimensions for the matrix
#'
#' @returns A downsampled matrix for rasterization
#'
#' @export

downsample_cor <- function(

  sub_mat,
  max_dim

  ) {
  n <- nrow(sub_mat)
  if (n <= max_dim) return(sub_mat)

  bin <- pmin(findInterval(seq_len(n), seq(1, n, length.out = max_dim + 1),
                           rightmost.closed = TRUE), max_dim)
  bin_size <- tabulate(bin, nbins = max_dim)

  row_avg <- rowsum(sub_mat, bin) / bin_size
  binned <- rowsum(t(row_avg), bin) / bin_size

  block_names <- sprintf("b%03d", seq_len(max_dim))
  dimnames(binned) <- list(block_names, block_names)
  binned
}


#' Prepare time data
#'
#' This function prepares time expression for group-mean plotting
#'
#' @param input_expression A matrix or data.frame object with rows and genes and columns as time points
#'
#' @returns A data.frame object formatted for use in [tiABBA::time_plot]
#'
#' @export

prepare_time <- function(

  input_expression

) {
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


#' Find IQR
#'
#' This function assists with quickly finding the IQR of a vector
#'
#' @param input_vector A numeric vector object
#'
#' @returns Two values containing the IQR - 1.5Q1 and IQR + 1.5Q3
#'
#' @export

find_iqr <- function(input_vector){
  iqr <- IQR(input_vector)
  q1 <- quantile(input_vector, 0.25)
  q3 <- quantile(input_vector, 0.75)

  low_lim <- (q1 - (1.5 * iqr))
  up_lim <- (q3 + (1.5 * iqr))

  return(c(low_lim, up_lim))
}
