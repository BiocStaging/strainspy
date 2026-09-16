testthat::test_that("strainspy_fit accessors return coherent values", {
  fit <- ss_fit_zib()

  testthat::expect_true(methods::validObject(fit))

  zi <- getZICoefficients(fit)
  testthat::expect_equal(nrow(zi), nrow(fit@coefficients))

  contigs <- getContigNames(fit)
  testthat::expect_length(contigs, nrow(fit@coefficients))

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
