#' @keywords internal
"_PACKAGE"

## R.utils is not called directly: data.table::fread() loads it to read the
## gzipped profile and taxonomy files. Importing one symbol keeps the
## dependency explicit and verified at install time.
## usethis namespace: start
#' @importFrom data.table :=
#' @importFrom data.table .BY
#' @importFrom data.table .EACHI
#' @importFrom data.table .GRP
#' @importFrom data.table .I
#' @importFrom data.table .N
#' @importFrom data.table .NGRP
#' @importFrom data.table .SD
#' @importFrom data.table data.table
#' @importFrom data.table set
#' @importFrom graphics abline
#' @importFrom methods slotNames
#' @importFrom stats binomial
#' @importFrom stats p.adjust
#' @importFrom stats prcomp
#' @importFrom stats rbinom
#' @importFrom utils capture.output
#' @importFrom ggrepel geom_text_repel
#' @importFrom harmonicmeanp p.hmp
#' @importFrom R.utils gunzip
#' @importFrom stringr str_detect
#' @importFrom tools file_path_sans_ext
## usethis namespace: end
NULL
