#' Prepare app data
#'
#' This function is a wrapper for other functions to create app data based off initial data
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param cor_seq A sequence of correlation filters to run on the data. Default value is the recommended range.
#' @param sil_seq A sequence of silhouette filters to run on the data. Default value is the recommended range.
#' @param k_max The maximum cluster amount to test. Results greater than 8 clusters may break plotting, along with being must more expensive to run.
#' @param parallel_silhouette FALSE or an integer setting the number of cores to use in parallel processing. Recommended to use [parallel::detectCores] and not exceed free memory divided by maximum object size.
#' @param set_workers An integer to set the amount of parallel workers used in plot generation.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns Locally saved results to data/group_condition and images to images/group_condition
#'
#' @export

wrap_app_inputs <- function(

  expr,
  group_file,
  cor_seq = seq(.1, .4, by = .01),
  sil_seq = seq(-.5, .7, by = .01),
  k_max = 8L,
  parallel_silhouette = FALSE,
  set_workers = {detectCores() - 1},
  verbose = TRUE

) {

  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  max_count <- nrow(expr)
  prepare_directory(group_file)

  cor_res <- scan_correlation(
    expr = expr,
    cor_seq = cor_seq,
    group_file = group_file,
    verbose = verbose
  )

  sil_res <- scan_silhouette(
    expr = expr,
    cor_res = cor_res,
    sil_seq = sil_seq,
    group_file = group_file,
    k_max = k_max,
    run_parallel = parallel_silhouette,
    verbose = verbose
  )

  scan_index(
    sil_res = sil_res,
    group_file = group_file,
    max_count = max_count
  )

  scan_plot(
    expr = expr,
    sil_res = sil_res,
    group_file = group_file,
    set_workers = set_workers,
    verbose = verbose
  )

  return(paste("Finished generating app data for", group_file))
}
