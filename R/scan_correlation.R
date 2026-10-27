#' Iterate Correlation Filters
#'
#' This function implements [totoInferno::cor_filter] across a parameter sequence. Saves intermediate results
#'
#' @param expr A matrix of expression data. Rows as genes and columns as time points.
#' @param cor_seq A sequence of correlation filters to run on the data. Default value is the recommended range.
#' @param group_file String for directory for the group_condition combination. Format as group_condition.
#' @param verbose Optional TRUE/FALSE to print progress along iterations.
#'
#' @returns A list object containing correlation filters tested and the genes kept. Duplicated results are removed.
#'
#' @export

scan_correlation <- function(

  expr,
  cor_seq = seq(.1, .4, by = .01),
  group_file,
  verbose = TRUE

) {
  if(!"matrix" %in% class(expr)){stop("Error: Please input a matrix class object")}

  expr <- as.matrix(expr)
  cor_res <- list()
  cor_data <- expr
  if(verbose == TRUE){
    cat("Running correlation scan: ")
    progress_bar <- txtProgressBar(min = 0, max = length(cor_seq), initial = 0, style = 3)
  }
  for(i in seq_along(cor_seq)){
    res_i <- tryCatch(cor_filter(expr = cor_data,
                                 cor_threshold = cor_seq[i],
                                 return_output = "expr"),

                      error = function(e) NULL)

    if(verbose == TRUE){setTxtProgressBar(progress_bar, i)}
    if(!is.null(res_i)){
      cor_res[[str_c(cor_seq[i])]] <- row.names(res_i)
      cor_data <- cor_data[row.names(res_i), , drop = FALSE]
    }else{break}
  }
  if(verbose == TRUE){close(progress_bar)}
  rm(cor_data)


  # filter to exclude duplicates
  cor_keep <- lengths(cor_res) %>% diff()
  cor_keep <- cor_keep[cor_keep != 0]

  cor_res <- cor_res[c("0.1", names(cor_keep))]

  saveRDS(cor_res, file.path("app/data", group_file, "cor_results"))
  return(cor_res)
}
