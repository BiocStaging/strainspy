testthat::test_that("glmZiBFit parallel backend matches the serial fit", {
  d <- ss_test_data()
  serial <- ss_fit_zib()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  parallel <- suppressWarnings(
    glmZiBFit(d$se_filtered[1:6, ], design, nthreads = 2, progress = FALSE)
  )

  testthat::expect_s4_class(parallel, "strainspy_fit")
  testthat::expect_equal(nrow(getCoefficients(parallel)), 6)

  # Preset priors depend only on the design, so per-row fits are independent
  # of which other rows are present: the parallel fit must reproduce the
  # cached serial fit exactly
  for (col in colnames(getCoefficients(parallel))) {
    testthat::expect_equal(getCoefficients(parallel)[[col]],
                           getCoefficients(serial)[[col]][1:6])
    testthat::expect_equal(getPValues(parallel)[[col]],
                           getPValues(serial)[[col]][1:6])
  }
})

testthat::test_that("glmZiBFit gamlss backend agrees with glmmTMB", {
  testthat::skip_if_not_installed("gamlss")
  testthat::skip_if_not_installed("gamlss.dist")

  s <- ss_sim()
  rows <- s$planted[1:4]

  fit_g <- suppressWarnings(
    glmZiBFit(s$se[rows, ], ~ Case_status, method = "gamlss",
              nthreads = 1, progress = FALSE)
  )

  testthat::expect_s4_class(fit_g, "strainspy_fit")
  testthat::expect_equal(nrow(getCoefficients(fit_g)), length(rows))
  testthat::expect_true(all(getConvergence(fit_g)))
  testthat::expect_true(all(getCoefficients(fit_g)$Case_statusPD < 0))
  testthat::expect_equal(getCoefficients(fit_g)$Case_statusPD,
                         getCoefficients(s$fit)$Case_statusPD[rows],
                         tolerance = 1e-2)
})

testthat::test_that("glmZiBFit rejects malformed inputs", {
  d <- ss_test_data()
  se <- d$se_filtered[1:3, ]

  testthat::expect_error(
    glmZiBFit(as.matrix(SummarizedExperiment::assay(se)), ~ Case_status,
              progress = FALSE),
    "must be a SummarizedExperiment"
  )
  testthat::expect_error(
    glmZiBFit(se, "~ Case_status", progress = FALSE),
    "must be a formula"
  )
  testthat::expect_error(
    glmZiBFit(se, ~ (1 | Case_status), progress = FALSE),
    "not a valid formula"
  )
  testthat::expect_error(
    glmZiBFit(se, ~ Case_status, method = "bogus", progress = FALSE),
    "should be one of"
  )
})
