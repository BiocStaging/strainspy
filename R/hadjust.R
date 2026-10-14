#' Performs hierarchical taxonomy aware adjustment of p-values.
#'
#' This function performs multiple testing adjustment on the results of a `strainspy_fit`
#' object at each taxonomic level, using either the Harmonic Mean P-value (HMP) or
#' Bonferroni correction.
#'
#' @param object A `strainspy_fit` object.
#' @param coef The coefficient to use for p-value adjustment. Defaults to 2.
#' @param method The method for p-value adjustment: "HMP" (default) or "bonferroni".
#' @param taxonomy Required. A taxonomy table from `read_taxonomy()`, or the path to the taxonomy file.
#' @param index_range Logical. Also return `index_min` and `index_max` columns, used to draw
#'   Manhattan plot segments for higher taxonomic levels. Default `FALSE`.
#' @return A tibble with one row per taxon and level, with columns `Level`, `Model`, `Name`,
#'   the number of strains (`strain_count` for HMP, `sequence_count` for Bonferroni),
#'   `mean_coefficient` and the adjusted p-value `p_adjust`.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' se <- filter_by_presence(se, min_nonzero = 30)
#' tax <- read_taxonomy(system.file("extdata", "example_taxonomy.tsv.gz",
#' package = "strainspy"))
#' fit <- glmZiBFit(se[1:10, ], ~ Case_status + Age_at_collection, 
#' nthreads = 1, return_vcov = FALSE)
#' out <- hadjust(fit, coef = 2, taxonomy = tax)
#' head(out)
#'
#' @importFrom tibble add_column as_tibble
#' @importFrom dplyr group_by left_join arrange mutate across all_of n rename summarise
#' @importFrom purrr map_dfr
#' @importFrom stats setNames
#'
#' @export
hadjust <- function(object, coef=2, method = "HMP", taxonomy, index_range=FALSE) {
  # Validate input
  if (!inherits(object, "strainspy_fit")) {
    stop("Input must be a strainspy_fit object.")
  }
  
  # Taxonomy is required: a table from read_taxonomy() or a path to the file
  if (is.null(taxonomy)) {
    stop("`taxonomy` is required.")
  } else if (is.character(taxonomy)) {
    tax <- read_taxonomy(taxonomy)
  } else {
    tax <- taxonomy
  }

  # Only two taxonomy-aware adjustments are implemented
  if (!method %in% c("HMP", "bonferroni")) {
    warning("method = '", method, "' is not available for taxonomy-aware adjustment; ",
            "using 'bonferroni'. Supported methods are 'HMP' and 'bonferroni'.", call. = FALSE)
  }
  
  # Extract p-values & coefficients
  mdl <- as.character(getFitCall(object))[[1]]
  if ("caseControlFit" == mdl) {
    main_model <- "Logistic"
  } else if ("glmZiBFit" == mdl | "glmObFit" == mdl) {
    main_model <- "Beta"
  } else {
    main_model <- "Unknown model" # We need to expand this as we include models
  }
  
  beta_res <- tibble::as_tibble(getRowData(object)) |>
    tibble::add_column(coefficient=getCoefficients(object)[[coef]]) |>
    tibble::add_column(std_error=getStdErrors(object)[[coef]]) |>
    tibble::add_column(p_value=getPValues(object)[[coef]]) |>
    tibble::add_column(Model=main_model, .before=1)
  
  if (!is.null(getZICoefficients(object))) {
    zi_res  <- tibble::as_tibble(getRowData(object)) |>
      tibble::add_column(coefficient=getZICoefficients(object)[[coef]]) |>
      tibble::add_column(std_error=getZIStdErrors(object)[[coef]]) |>
      tibble::add_column(p_value=getZIPValues(object)[[coef]]) |>
      tibble::add_column(Model="Zero-Inflated", .before=1)
  }
  
  # Check which column in row_data is corresponds to the taxonomic data
  row_data <- getRowData(object)
  tax_col <- Position(function(i) all(row_data[, i] %in% tax$Genome),
                      seq_len(ncol(row_data)))
  if (is.na(tax_col)) {
    stop("Taxon names in `taxonomy` do not match any column of rowData(object); ",
         "cannot adjust p-values by taxonomy.")
  }
  
  hierarchy_df <- beta_res |>
    dplyr::left_join(tax, by = setNames("Genome", names(getRowData(object))[tax_col])) |>
    dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index) |>
    dplyr::mutate(index=seq_len(nrow(beta_res)))
  
  if (!is.null(getZICoefficients(object))){
    hierarchy_df <- rbind(
      hierarchy_df,
      zi_res |>
        dplyr::left_join(tax, by = setNames("Genome", names(getRowData(object))[tax_col]))  |>
        dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index) |>
        dplyr::mutate(index=seq_len(nrow(zi_res)))
    )
  }
  
  # Remove rows which did not converge or errored
  hierarchy_df <- hierarchy_df[!is.na(hierarchy_df$coefficient),]
  
  if(!is.null(tax_col)){ 
    # Now reorder by global taxonomic tree
    ordered_genomes <- colour_by_tax(unique(hierarchy_df$Genome_file), tax)
    
    hierarchy_df <- hierarchy_df %>%
      dplyr::arrange(Model, match(Genome_file, ordered_genomes)) %>%
      dplyr::group_by(Model) %>%
      dplyr::mutate(index = dplyr::dense_rank(match(Genome_file, ordered_genomes))) %>%
      dplyr::ungroup()
  }
  # Define the grouping columns
  grouping_cols <- c("Phylum", "Class", "Order", "Family", "Genus", "Species", colnames(hierarchy_df)[[3]])
  
  if (method == "HMP") {
    L <- nrow(beta_res)
    hierarchy_df$w <- 1/nrow(beta_res)
    
    
    adj <- purrr::map_dfr(grouping_cols, ~{
      hierarchy_df |>
        dplyr::group_by(Model, dplyr::across(dplyr::all_of(.x))) |>  # Use all_of to refer to the grouping column
        dplyr::summarise(
          strain_count = dplyr::n(),
          mean_coefficient = mean(coefficient),
          p_adjust = p.hmp(p_value, w, L, multilevel = TRUE) / sum(w),
          index_min = min(index),
          index_max = max(index),
          .groups = "drop"  # Drop the grouping after summarizing
        ) |>
        dplyr::rename(Name = all_of(.x)) |>  # Use all_of to rename the grouped column
        tibble::add_column(Level = .x, .before = 1)  # Add the 'Level' column before the first column
    }) |>
      dplyr::arrange(p_adjust)
  } else {
    hierarchy_df <- hierarchy_df |>
      dplyr::group_by(Model) |>
      dplyr::mutate(padj = p.adjust(p_value, method = "bonferroni")) |>
      dplyr::ungroup()
    
    adj <- purrr::map_dfr(grouping_cols, ~{
      hierarchy_df |>
        dplyr::group_by(Model, dplyr::across(dplyr::all_of(.x))) |>  # Use all_of to refer to the grouping column
        dplyr::summarise(
          sequence_count = dplyr::n(),
          mean_coefficient = mean(coefficient),
          p_adjust = min(padj),
          index_min = min(index),
          index_max = max(index),
          .groups = "drop"  # Drop the grouping after summarizing
        ) |>
        dplyr::rename(Name = all_of(.x)) |>  # Use all_of to rename the grouped column
        tibble::add_column(Level = .x, .before = 1)  # Add the 'Level' column before the first column
    }) |>
      dplyr::arrange(p_adjust)
  }
  
  if(!index_range) {
    adj$index_min <- NULL
    adj$index_max <- NULL
  }
  
  return(adj)
}
