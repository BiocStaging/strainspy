testthat::test_that("strainspy_fit accessors return coherent values", {
  fit <- ss_fit_zib()

  testthat::expect_true(methods::validObject(fit))

  zi <- getZICoefficients(fit)
  testthat::expect_equal(nrow(zi), nrow(getCoefficients(fit)))

  contigs <- getContigNames(fit)
  testthat::expect_length(contigs, nrow(getCoefficients(fit)))

  testthat::expect_true(length(getGenomes(fit)) > 0)
  testthat::expect_true(length(getSamples(fit)) > 0)

  res <- getResiduals(fit)
  testthat::expect_false(is.null(res))
})

testthat::test_that("strainspy_fit show and summary methods run", {
  fit <- ss_fit_zib()

  testthat::expect_output(show(fit))
  testthat::expect_output(print(fit))

  s <- summary(fit)
  testthat::expect_false(is.null(s))
})

testthat::test_that("every strainspy_fit accessor returns the documented type", {
  fit_z <- ss_fit_zib()
  n <- nrow(getRowData(fit_z))

  for (acc in c(getCoefficients, getStdErrors, getPValues, getZIStdErrors, getZIPValues)) {
    testthat::expect_s4_class(acc(fit_z), "DFrame")
    testthat::expect_equal(nrow(acc(fit_z)), n)
  }
  testthat::expect_s4_class(getVcov(fit_z), "List")
  testthat::expect_s4_class(getZIVcov(fit_z), "List")
  testthat::expect_s3_class(getDesign(fit_z), "formula")
  testthat::expect_s3_class(getFamily(fit_z), "family")
  testthat::expect_type(getFitCall(fit_z), "language")
  testthat::expect_true(is.logical(isScaleContinuous(fit_z)))
  testthat::expect_length(getConvergence(fit_z), n)
})
