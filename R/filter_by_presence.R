#' filter_by_presence
#'
#' This function filters a SummarizedExperiment object to keep only rows (features/strains)
#' that are present in a specified number of samples.
#' The same rows are filtered from all assays in the SummarizedExperiment.
#'
#' @param se A `SummarizedExperiment` object.
#' @param min_nonzero Integer. The minimum number of non-zero entries required to retain a row. Default is 10.
#'
#' @return A filtered `SummarizedExperiment` object with only the rows that meet the criteria.
#'
#' @import SummarizedExperiment
#'
#' @examples
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy")
#' se <- read_sylph(example_path)
#'
#' # Filter to keep only rows with at least 20 non-zero entries 
#' filtered_se <- filter_by_presence(se, min_nonzero = 20)
#'
#' @export
filter_by_presence <- function(se, min_nonzero = 10) {

  # Check that the input is a SummarizedExperiment object
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment object.")
  }

  if(min_nonzero %% 1 != 0) {
    stop("`min_nonzero` must be an integer.")
  }

  if(min_nonzero > dim(se)[2]){
    stop("Ensure `min_nonzero` <= number of samples in `se`: ", dim(se)[2])
  }
  # Count the number of non-zero entries in each row of the assay
  nonzero_counts <- Matrix::rowSums(SummarizedExperiment::assays(se)[[1]] != 0)
  
  # Identify which rows have at least 'min_nonzero' non-zero entries
  rows_to_keep <- nonzero_counts >= min_nonzero
    
  # Filter the assays, rowData, and colData in the SummarizedExperiment
  message("Retained ", sum(rows_to_keep), " rows after filtering")
  filtered_se <- se[rows_to_keep, ]
  
  # Return the filtered SummarizedExperiment
  return(filtered_se)
}
