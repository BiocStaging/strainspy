ss_test_data <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    meta_path <- system.file("extdata", "example_metadata.csv.gz", package = "strainspy")
    sylph_path <- system.file("extdata", "example_sylph_profile.tsv.gz", package = "strainspy")
    taxonomy_path <- system.file("extdata", "example_taxonomy.tsv.gz", package = "strainspy")
    metaphlan_path <- system.file("extdata", "metaphlan_merged.tsv.gz", package = "strainspy")
    sourmash_path <- system.file("extdata", "example_sourmash.csv.gz", package = "strainspy")

    stopifnot(
      nzchar(meta_path),
      nzchar(sylph_path),
      nzchar(taxonomy_path),
      nzchar(metaphlan_path),
      nzchar(sourmash_path)
    )

    meta <- read.csv(meta_path, stringsAsFactors = FALSE)
    meta$Case_status <- as.factor(meta$Case_status)

    se <- strainspy::read_sylph(sylph_path, meta_data = meta)
    se_filtered <- strainspy::filter_by_presence(se, min_nonzero = 30)
    taxonomy <- strainspy::read_taxonomy(taxonomy_path)

    cache <<- list(
      meta = meta,
      meta_path = meta_path,
      sylph_path = sylph_path,
      taxonomy_path = taxonomy_path,
      metaphlan_path = metaphlan_path,
      sourmash_path = sourmash_path,
      se = se,
      se_filtered = se_filtered,
      taxonomy = taxonomy
    )

    cache
  }
})

ss_fit_zib <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- suppressWarnings(
      strainspy::glmZiBFit(
        d$se_filtered[1:20, ],
        design,
        nthreads = 1,
        return_vcov = TRUE
      )
    )

    cache
  }
})

ss_fit_zib_readme <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- suppressWarnings(
      strainspy::glmZiBFit(
        d$se_filtered,
        design,
        nthreads = 1,
        return_vcov = FALSE
      )
    )

    cache
  }
})

ss_fit_qb <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- strainspy::glmQBFit(
      d$se_filtered[1:30, ],
      design,
      nthreads = 1
    )

    cache
  }
})

ss_priors_eb <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    design <- stats::as.formula("~ Case_status + Age_at_collection")

    cache <<- suppressWarnings(
      strainspy::compute_eb_priors(
        d$se_filtered[1:10, ],
        design,
        nthreads = 1,
        progress = FALSE
      )
    )

    cache
  }
})

# Simulated ANI data with a known planted effect. Rows `planted` have lower
# identity and more zeros in PD than in Control; rows `nulls` are drawn from
# the same distribution in both groups. Built on top of a real fixture so the
# SummarizedExperiment carries the rowData/colData the fitters expect.
ss_sim <- local({
  cache <- NULL

  function() {
    if (!is.null(cache)) {
      return(cache)
    }

    d <- ss_test_data()
    se <- d$se_filtered[1:12, ]
    grp <- SummarizedExperiment::colData(se)$Case_status
    planted <- 1:6
    nulls <- 7:12

    # ~98% identity, present in roughly half the samples (values in [0, 1])
    set.seed(42)
    base <- function(n) stats::rbeta(n, 200, 4) * stats::rbinom(n, 1, 0.5)

    mat <- matrix(0, nrow(se), ncol(se),
                  dimnames = dimnames(SummarizedExperiment::assay(se)))
    for (i in planted) {
      mat[i, grp == "Control"] <- base(sum(grp == "Control")) * 100
      mat[i, grp == "PD"] <- strainspy:::rescale_beta(
        base(sum(grp == "PD")), beta = 0.99, zi = 0.15
      )$rescaled * 100
    }
    for (i in nulls) {
      mat[i, ] <- base(ncol(se)) * 100
    }
    SummarizedExperiment::assay(se, 1) <- mat

    fit <- suppressWarnings(
      strainspy::glmZiBFit(se, ~ Case_status, nthreads = 1, progress = FALSE)
    )

    cache <<- list(se = se, fit = fit, planted = planted, nulls = nulls)
    cache
  }
})
