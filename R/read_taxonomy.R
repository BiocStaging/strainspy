#' Read and process a taxonomy file
#'
#' This function reads a headerless, two-column, tab-separated taxonomy file: genome
#' accession and GTDB-style lineage string (`d__...;p__...;...;s__...`), as distributed
#' with the sylph GTDB databases.
#'
#' @importFrom data.table fread tstrsplit
#' @importFrom dplyr arrange
#' 
#' @param file_path Path to the taxonomy file.
#' @return A data.table with columns `Genome`, `Taxonomy`, `Domain`, `Phylum`, `Class`,
#'   `Order`, `Family`, `Genus`, `Species` and `index` (the row order after sorting by lineage).
#' @examples
#' tax <- read_taxonomy(system.file("extdata", "example_taxonomy.tsv.gz", 
#' package = "strainspy"))
#' head(tax)
#' @export
read_taxonomy <- function(file_path) {
  # Read the taxonomy file
  tax <- data.table::fread(file_path, sep = "\t", col.names = c("Genome", "Taxonomy"))

  # Clean and split the Taxonomy column into individual taxonomic ranks
  tax$Taxonomy <- gsub("[a-z]__", "", tax$Taxonomy)
  tax[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") := data.table::tstrsplit(Taxonomy, ";", fill = TRUE)]

  # Add index and arrange
  tax$index <- seq_len(nrow(tax))
  tax <- tax |>
    dplyr::arrange(Domain, Phylum, Class, Order, Family, Genus, Species, index)
  tax$index <- seq_len(nrow(tax))

  return(tax)
}

