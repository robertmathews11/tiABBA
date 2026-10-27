#' Iterate Silhouette Filters
#'
#' This function implements [totoInferno::silhouette_filter] across a parameter sequence. Saves intermediate results.
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param cor_res A list containing the results from [tiABBA::scan_correlation]
#' @param sil_seq A sequence of silhouette filters to run on the data. Default value is the recommended range.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param k_max The maximum cluster amount to test. Results greater than 8 clusters may break plotting, along with being must more expensive to run.
#' @param run_parallel FALSE or an integer setting the number of cores to use in parallel processing. Recommended to use [parallel::detectCores] and not exceed free memory divided by maximum object size.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns A list object with multiple levels: correlation parameter, silhouette parameter, then results
#'
#' @export

scan_silhouette <- function(

  expr,
  cor_res,
  sil_seq = seq(-.5, .7, by = .01),
  group_file,
  k_max = 8L,
  run_parallel = FALSE,
  verbose = TRUE

) {
  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  ## Run silhouette, output panels
  sil_seq = seq(-.5, .7, by = .01)

  if(verbose == TRUE){cat("Running silhouette scan: ")}
  if(run_parallel == FALSE){
    sil_res <- imap(cor_res, \(cor_i, cor_name) {
      cor_filter <- expr[cor_i, , drop = FALSE]
      pre_cor <- cor(t(cor_filter), method = "pearson")

      sil_res <- list()
      for(i in seq_along(sil_seq)){
        res_i <- tryCatch(silhouette_filter(expr = cor_filter,
                                            sil_threshold = sil_seq[i],
                                            cor_expr = pre_cor,
                                            k_max = k_max),

                          error = function(e) NULL)

        if(!is.null(res_i)){

          sil_res[[str_c(sil_seq[i])]] <- filter_res(res_i)

          cor_filter <- cor_filter[row.names(res_i[["expr"]]), , drop = FALSE]
        }else{break}
      }
      return(sil_res)
    }, .progress = {verbose})
  }else{
    plan(multisession, workers = run_parallel)
    sil_res <- future_imap(cor_res, \(cor_i, cor_name) {
      cor_filter <- expr[cor_i, , drop = FALSE]
      pre_cor <- cor(t(cor_filter), method = "pearson")

      sil_res <- list()
      for(i in seq_along(sil_seq)){
        res_i <- tryCatch(silhouette_filter(expr = cor_filter,
                                            sil_threshold = sil_seq[i],
                                            cor_expr = pre_cor,
                                            k_max = 8L),

                          error = function(e) NULL)

        if(!is.null(res_i)){

          sil_res[[str_c(sil_seq[i])]] <- filter_res(res_i)

          cor_filter <- cor_filter[row.names(res_i[["expr"]]), , drop = FALSE]
        }else{break}
      }
      return(sil_res)
    }, .progress = {verbose})
    plan(sequential)
  }

  sil_res <- map(sil_res, \(res_i){
    correct_names <- names(res_i) %>%
      as.numeric() %>%
      round(digits = 2)
    return(setNames(res_i, correct_names))
  })

  saveRDS(sil_res, file = file.path("app/data", group_file, "results"))
  return(sil_res)
}
