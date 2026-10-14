#' @importFrom S4Vectors DataFrame List

setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))
setClassUnion("strainspyPriorsOrNULL", c("strainspy_priors", "NULL"))
setClassUnion("ListOrNULL", c("List", "NULL"))

#' strainspy_fit: A Class for storing fit data from all strainspy models
#'
#' The `strainspy_fit` class stores the results of a strainspy model fit.
#' 
#' @importClassesFrom S4Vectors DataFrame List
#' @importMethodsFrom S4Vectors show
#'
#'
#' @slot row_data A `DataFrame` containing row data.
#' @slot priors A `strainspy_priors` object, or `NULL` for an unpenalised maximum-likelihood fit.
#' @slot scale_continuous A logical indicating whether all numeric columns were z-score standardized
#' @slot coefficients A `DataFrame` containing model coefficients.
#' @slot std_errors A `DataFrame` containing standard errors of the coefficients.
#' @slot p_values A `DataFrame` containing p-values for the coefficients.
#' @slot vcov A `List` of per-feature variance-covariance matrices, or `NULL`.
#' @slot zi_coefficients A `DataFrame` containing zero-inflation coefficients.
#' @slot zi_std_errors A `DataFrame` containing zero-inflation standard errors.
#' @slot zi_p_values A `DataFrame` containing zero-inflation p-values.
#' @slot zi_vcov A `List` of per-feature zero-inflation variance-covariance matrices, or `NULL`.
#' @slot residuals A `DataFrame` containing residuals from the model.
#' @slot convergence A logical indicating whether model fitting converged.
#' @slot design The model formula.
#' @slot call The matched call of the model.
#' @slot family The family of the fit. 
#'
#' @export
methods::setClass("strainspy_fit",
                  slots = list(
                    row_data = "DataFrame",
                    priors = "strainspyPriorsOrNULL",
                    scale_continuous = "logical",
                    coefficients = "DataFrame",
                    std_errors = "DataFrame",
                    p_values = "DataFrame",
                    vcov = "ListOrNULL",
                    zi_coefficients = "DataFrameOrNULL",
                    zi_std_errors = "DataFrameOrNULL",
                    zi_p_values = "DataFrameOrNULL",
                    zi_vcov = "ListOrNULL",
                    residuals = "DataFrameOrNULL",
                    convergence = "logical",
                    design = "formula",
                    call = "call",
                    family = "ANY" # We should make this mandatory
                  )
)

# Every per-feature slot is indexed positionally against `row_data`, so a length
# mismatch silently attributes a strain's name to another strain's statistics.
# Fail loudly instead.
methods::setValidity("strainspy_fit", function(object) {
  n <- nrow(object@row_data)
  per_feature <- list(
    coefficients = object@coefficients,
    std_errors = object@std_errors,
    p_values = object@p_values,
    zi_coefficients = object@zi_coefficients,
    zi_std_errors = object@zi_std_errors,
    zi_p_values = object@zi_p_values,
    vcov = object@vcov,
    zi_vcov = object@zi_vcov
  )
  bad <- vapply(per_feature, function(x) {
    if (is.null(x)) return(FALSE)
    len <- if (is(x, "DataFrame")) nrow(x) else length(x)
    len != n
  }, logical(1))
  if (length(object@convergence) && length(object@convergence) != n) {
    bad <- c(bad, convergence = TRUE)
  }
  if (any(bad)) {
    return(paste0(
      "slot(s) ", paste(names(bad)[bad], collapse = ", "),
      " are not aligned with row_data (", n, " features)."
    ))
  }
  TRUE
})

#' Show Method for strainspy_fit
#'
#' Displays a summary of a `strainspy_fit` object, including model design,
#' convergence, and coefficient matrix dimensions.
#'
#' @param object An object of class `strainspy_fit`.
#' @return Invisibly returns `NULL` after printing a compact summary.
#' @export
methods::setMethod("show", "strainspy_fit", function(object) {
  cat("strainspy_fit object\n")
  cat(rep("-", 30), "\n", sep = "")
  
  cat("Design formula: ")
  print(object@design)
  
  cat(rep("-", 30), "\n", sep = "")
  cat("Converged=", length(which(object@convergence==TRUE)) , ", Failed=", length(which(object@convergence==FALSE)), "\n", sep = '')
  
  cat(rep("-", 30), "\n", sep = "")
  cat("Coefficients: ", nrow(object@coefficients), " features x ", ncol(object@coefficients), " terms\n")
  
  cat(rep("-", 30), "\n", sep = "")
  if (!is.null(object@zi_coefficients)) {
    cat("ZI Coefficients: ", nrow(object@zi_coefficients), " features x ", ncol(object@zi_coefficients), " terms\n")
  } else {
    cat("ZI Coefficients: None\n")
  }
  
  cat(rep("-", 30), "\n", sep = "")
  cat("Residuals available: ", ifelse(!is.null(object@residuals), "Yes", "No"), "\n")
  cat("Call:\n")
  print(object@call)
  cat(rep("-", 30), "\n", sep = "")
  
  if(!is.null(object@priors)) {
    cat("Priors used:\n")
    show(object@priors)
  } else {
    cat("Fitted without MAP priors (maximum likelihood).\n")
  }
  
  cat(rep("-", 30), "\n", sep = "")
})

#' Access Zero-Inflation Coefficients
#'
#' Returns the zero-inflation coefficients from a `strainspy_fit` object.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A `DataFrame` containing zero-inflation coefficients.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' class(getZICoefficients(fit))
#' @export
methods::setGeneric("getZICoefficients", function(object) methods::standardGeneric("getZICoefficients"))

#' @rdname getZICoefficients
#' @export
methods::setMethod("getZICoefficients", "strainspy_fit", function(object) {
  object@zi_coefficients
})

#' Access Residuals
#'
#' Returns the residuals from a `strainspy_fit` object.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A `DataFrame` containing residuals.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' class(getResiduals(fit))
#' @export
methods::setGeneric("getResiduals", function(object) methods::standardGeneric("getResiduals"))

#' @rdname getResiduals
#' @export
methods::setMethod("getResiduals", "strainspy_fit", function(object) {
  object@residuals
})

#' Summary Method for strainspy_fit
#'
#' Provides a high-level summary of a `strainspy_fit` object: its call, the first
#' rows of the coefficients and zero-inflation coefficients, and whether residuals
#' are available.
#'
#' @param object An object of class `strainspy_fit`.
#' @return A list summarizing the `strainspy_fit` object.
#'
#' @importFrom utils head
#'
#' @export
methods::setMethod("summary", "strainspy_fit", function(object) {
  list(
    call = object@call,
    coefficients = head(object@coefficients),
    zi_coefficients = head(object@zi_coefficients),
    residuals_available = !is.null(object@residuals)
  )
})


#' Get Contig names in a fit object
#'
#' Gets the ordered list of fitted contig names in a `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of contig names
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' head(getContigNames(fit))
#' @export
methods::setGeneric("getContigNames", function(object) methods::standardGeneric("getContigNames"))

#' @rdname getContigNames
#' @export
methods::setMethod("getContigNames", "strainspy_fit", function(object) {
  object@row_data$Contig_name
})

#' Get Genomes in a fit object
#'
#' Gets the ordered list of fitted genomes in a `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of fitted genomes
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' head(getGenomes(fit))
#' @export
methods::setGeneric("getGenomes", function(object) methods::standardGeneric("getGenomes"))

#' @rdname getGenomes
#' @export
methods::setMethod("getGenomes", "strainspy_fit", function(object) {
  object@row_data$Genome_file
})

#' Get Sample names in a fit object
#'
#' Gets the ordered list of samples used in `strainspy_fit`
#'
#' @param object An object of class `strainspy_fit`.
#' @return A vector of fitted sample names 
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ], 
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' head(getSamples(fit))
#' @export
methods::setGeneric("getSamples", function(object) methods::standardGeneric("getSamples"))

#' @rdname getSamples
#' @export
methods::setMethod("getSamples", "strainspy_fit", function(object) {
  names(object@residuals)
})

#' Accessors for strainspy_fit objects
#'
#' Getters for the per-feature statistics and fitting metadata stored in a
#' `strainspy_fit` object. Use these rather than reaching into the object
#' with `@`.
#'
#' @param object An object of class `strainspy_fit`.
#' @return
#' * `getRowData()`: a `DataFrame` of feature annotations, one row per feature.
#' * `getCoefficients()`, `getStdErrors()`, `getPValues()`: `DataFrame`s of
#'   per-feature estimates, one column per model term.
#' * `getZIStdErrors()`, `getZIPValues()`: the zero-inflation counterparts, or
#'   `NULL` when the model has no zero-inflation component.
#' * `getVcov()`, `getZIVcov()`: a `List` of per-feature variance-covariance
#'   matrices, or `NULL` when they were not stored.
#' * `getDesign()`: the model formula.
#' * `getFamily()`: the model family object.
#' * `getFitCall()`: the call that produced the fit.
#' * `isScaleContinuous()`: `TRUE` if continuous covariates were z-scored
#'   before fitting.
#' * `getConvergence()`: a logical vector, one entry per feature, `TRUE` where
#'   the model converged.
#' @examples
#' meta <- read.csv(system.file("extdata", "example_metadata.csv.gz",
#' package = "strainspy"))
#' meta$Case_status <- factor(meta$Case_status)
#' se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz",
#' package = "strainspy"), meta_data = meta)
#' fit <- glmZiBFit(filter_by_presence(se, min_nonzero = 30)[1:10, ],
#' ~ Case_status + Age_at_collection, nthreads = 1, return_vcov = FALSE)
#' head(getCoefficients(fit))
#' getDesign(fit)
#' isScaleContinuous(fit)
#' @name strainspy_fit-accessors
#' @rdname strainspy_fit-accessors
NULL

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getRowData", function(object) methods::standardGeneric("getRowData"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getRowData", "strainspy_fit", function(object) object@row_data)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getCoefficients", function(object) methods::standardGeneric("getCoefficients"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getCoefficients", "strainspy_fit", function(object) object@coefficients)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getStdErrors", function(object) methods::standardGeneric("getStdErrors"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getStdErrors", "strainspy_fit", function(object) object@std_errors)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getPValues", function(object) methods::standardGeneric("getPValues"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getPValues", "strainspy_fit", function(object) object@p_values)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getZIStdErrors", function(object) methods::standardGeneric("getZIStdErrors"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getZIStdErrors", "strainspy_fit", function(object) object@zi_std_errors)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getZIPValues", function(object) methods::standardGeneric("getZIPValues"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getZIPValues", "strainspy_fit", function(object) object@zi_p_values)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getVcov", function(object) methods::standardGeneric("getVcov"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getVcov", "strainspy_fit", function(object) object@vcov)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getZIVcov", function(object) methods::standardGeneric("getZIVcov"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getZIVcov", "strainspy_fit", function(object) object@zi_vcov)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getDesign", function(object) methods::standardGeneric("getDesign"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getDesign", "strainspy_fit", function(object) object@design)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getFamily", function(object) methods::standardGeneric("getFamily"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getFamily", "strainspy_fit", function(object) object@family)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getFitCall", function(object) methods::standardGeneric("getFitCall"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getFitCall", "strainspy_fit", function(object) object@call)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("isScaleContinuous", function(object) methods::standardGeneric("isScaleContinuous"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("isScaleContinuous", "strainspy_fit", function(object) object@scale_continuous)

#' @rdname strainspy_fit-accessors
#' @export
methods::setGeneric("getConvergence", function(object) methods::standardGeneric("getConvergence"))
#' @rdname strainspy_fit-accessors
#' @export
methods::setMethod("getConvergence", "strainspy_fit", function(object) object@convergence)
