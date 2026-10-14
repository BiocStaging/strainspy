#' Extract data from a StrainSpy `SummarizedExperiment` object in long format
#' @importFrom SummarizedExperiment assay
#' @importFrom tidyr pivot_longer
#'
#' @param se  SummarizedExperiment. A `SummarizedExperiment` object containing the assay data and metadata.
#' @param variables One or multiple variables in `colnames(colData(se))`, except `Sample_file`.
#' @param contigs A character vector of contigs (rownames of `se`) to extract.
#' @param taxonomy An optional taxonomy file read using strainspy::read_taxonomy (default=NULL).
#' @param drop_zeros Bool. If TRUE, rows with an ANI or abundance of 0 are dropped. Default FALSE.
#' @return A long-format `data.frame` with sample metadata, contig identity columns,
#' ANI values, row-level annotations, and optional taxonomy columns.
#'
#' @examples
#' library(strainspy)
#'
#' example_meta_path <- system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy")
#' example_meta <- read.csv(example_meta_path)
#' # Change variable from character to factor
#' example_meta$Case_status = factor(example_meta$Case_status)
#'
#' example_path <- system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy")
#' se <- read_sylph(example_path, example_meta)
#'
#' tax_path <- system.file("extdata", "example_taxonomy.tsv.gz", 
#' package = "strainspy")
#' tax <- read_taxonomy(tax_path)
#'
#' extract_strainspy(se, variables = c("Case_status", "Age_at_collection"), 
#' contigs = rownames(se)[1:5], taxonomy=tax)
#'
#'
#' @export
extract_strainspy <- function(se, variables=NULL, contigs=NULL, taxonomy=NULL, drop_zeros = FALSE){
  
  # default variables to colnames of colData
  if (is.null(variables)) {
    variables <- colnames(SummarizedExperiment::colData(se))
  }
  
  if(! all(variables %in% colnames(SummarizedExperiment::colData(se))) ){
    stop("At least one variable not found in colnames(colData(se))")
  }
  
  # check contigs
  chk_contigs <- contigs %in% rownames(se)
  if( ! all(chk_contigs)  ){
    # Name the offending contigs in the condition itself. Printing them to
    # stdout put them on a different stream from the error on stderr, so a
    # caller could see "Above contigs" with nothing above it.
    missing_txt <- paste(contigs[!chk_contigs], collapse = ", ")
    stop("Contigs missing from rownames(se): ", missing_txt)
  }
  
  tmp <- tibble::rownames_to_column(
    as.data.frame(t(as.matrix(SummarizedExperiment::assay(se[contigs, ])))),
    var = "sample")
  tmp <- dplyr::bind_cols(as_tibble(SummarizedExperiment::colData(se)[variables]), tmp)
  
  df_long <- tidyr::pivot_longer(tmp, cols = -c(sample, dplyr::all_of(variables)), names_to = "Contig_name", values_to = "ANI")
  
  if (drop_zeros) {
    df_long <- df_long |> filter(ANI > 0)
  }
  
  # Add other columns from rowdata
  df_long <- df_long |>
    left_join(tibble::as_tibble(SummarizedExperiment::rowData(se)), by = "Contig_name")
  
  if (!is.null(taxonomy)) {
    if (!is.data.frame(taxonomy)) {
      stop("`taxonomy` must be a data.frame, e.g. the output of read_taxonomy().")
    }
    if (!"Genome" %in% names(taxonomy)) {
      stop("Column 'Genome' not found in `taxonomy`.")
    }
    if (!"Genome_file" %in% names(df_long)) {
      stop("Column 'Genome_file' not found in rowData(se); cannot join taxonomy.")
    }
    df_long <- df_long |>
      left_join(taxonomy, by = c("Genome_file" = "Genome"))
  }
  
  # Reorder columns in df_long
  df_long <- df_long |>
    dplyr::select(sample, dplyr::all_of(variables), Contig_name, dplyr::everything(), ANI)
  
  return(df_long)
}
