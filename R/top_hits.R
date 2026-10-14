#' Provides the top results of a strainspy_fit object, sorted by adjusted p-value.
#'
#' Adjusts the p-values of one coefficient of a `strainspy_fit` for multiple testing with
#' `p.adjust()` and returns the features that pass `alpha`, sorted by adjusted p-value.
#' For taxonomy-aware adjustment see `hadjust()`.
#'
#' @param fit A `strainspy_fit` object.
#' @param coef The coefficient to report and use for p-value adjustment. Defaults to 2.
#' @param method The method for p-value adjustment, passed to `p.adjust()`. One of
#'   "holm" (default), "bonferroni", "BH" (Benjamini-Hochberg) or "BY".
#' @param alpha Numeric. Significance threshold for adjusted p-values. Defaults to 0.05.
#' 
#' @return A tibble with the top hits, sorted by adjusted p-value.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz",
#'  package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' top_hits(fit, coef = 2, alpha = 1)
#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble add_column
top_hits <- function(fit, coef=2, method = "holm", alpha=0.05) {
  # Validate input
  if (!inherits(fit, "strainspy_fit")) {
    stop("Input must be a strainspy_fit")
  }

  if (!is.numeric(coef) || length(coef) != 1 || is.na(coef) || coef %% 1 != 0) {
    stop("`coef` must be a single integer index.")
  }

  if (coef < 1 || coef > ncol(getCoefficients(fit))) {
    stop(sprintf("`coef` is out of bounds. Choose a value between 1 and %d.", ncol(getCoefficients(fit))))
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single numeric value between 0 and 1.")
  }
  
  # Check method is one of the available options
  if (!method %in% c("bonferroni", "BH", "BY", "holm")) {
    stop("Method must be one of 'bonferroni', 'BH', 'BY', or 'holm'.")
  }
  
 
  
  # Extract p-values & coefficients
  res <- tibble::as_tibble(getRowData(fit)) |>
    tibble::add_column(coefficient=getCoefficients(fit)[[coef]]) |>
    tibble::add_column(std_error=getStdErrors(fit)[[coef]]) |>
    tibble::add_column(p_value=getPValues(fit)[[coef]]) |>
    tibble::add_column(p_adjust=p.adjust(getPValues(fit)[[coef]], method = method))
  
  # Add zero-inflated coefficients if available
  if (!is.null(getZICoefficients(fit))){
    res <- res |>
      tibble::add_column(zi_coefficient=getZICoefficients(fit)[[coef]]) |>
      tibble::add_column(zi_std_error=getZIStdErrors(fit)[[coef]]) |>
      tibble::add_column(zi_p_value=getZIPValues(fit)[[coef]]) |>
      tibble::add_column(zi_p_adjust=p.adjust(getZIPValues(fit)[[coef]], method = method))
  }
  
  ## Arrange by adjusted p_value
  if(!is.null(getZICoefficients(fit))) {
    res <- res |> dplyr::arrange(pmin(p_adjust, zi_p_adjust))
  } else {
    res <- res |> dplyr::arrange(p_adjust)
  }
  
  # Filter on provided alpha
  if (!is.null(getZICoefficients(fit))){
    res <- res |> dplyr::filter(pmin(p_adjust, zi_p_adjust) <= alpha)
  } else {
    res <- res |> dplyr::filter(p_adjust <= alpha)
  }
  
  if(nrow(res) == 0) {
    warning(sprintf("Multiple testing correction using `%s`: No significant associations detected for coef = %d at alpha = %f", method, coef, alpha))
  } else {
    message("Found ", nrow(res), " tophits for ", names(getCoefficients(fit))[coef],
            " at alpha = ", alpha, " using ", method)
    # shortcut to attach phenotype and other info to top_hits
    attr(res, "phenotype_coef") <- coef
    attr(res, "method") <- method
    attr(res, "alpha") <- alpha
  }
  
  
  
  return(res)
}
