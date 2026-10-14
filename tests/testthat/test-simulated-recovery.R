testthat::test_that("glmZiBFit recovers a planted ANI effect", {
  s <- ss_sim()
  fit <- s$fit

  testthat::expect_true(all(getConvergence(fit)))

  beta <- getCoefficients(fit)$Case_statusPD
  beta_p <- getPValues(fit)$Case_statusPD

  # Lower identity in PD must come back as a clearly negative coefficient
  testthat::expect_true(all(beta[s$planted] < 0))
  testthat::expect_true(all(beta_p[s$planted] < 0.01))

  # Null rows carry no signal
  testthat::expect_true(all(beta_p[s$nulls] > 0.05))

  # Every planted effect is more negative than every null estimate
  testthat::expect_lt(max(beta[s$planted]), min(beta[s$nulls]))

  # Extra zeros in PD shift the zero-inflation component upwards on average
  zi <- getZICoefficients(fit)$Case_statusPD
  testthat::expect_gt(mean(zi[s$planted]), mean(zi[s$nulls]))
})

testthat::test_that("top_hits detects the planted rows", {
  s <- ss_sim()
  contigs <- SummarizedExperiment::rowData(s$se)$Contig_name

  th <- top_hits(s$fit, coef = 2)
  cm <- strainspy:::get_confusion_mx(th, contigs[s$planted], contigs,
                                     print_cm = FALSE)

  # All planted rows recovered; at most one chance false positive, which can
  # arise from the zero-inflation component on a null row
  testthat::expect_equal(cm["Positive", "Positive"], 6)
  testthat::expect_equal(cm["Negative", "Positive"], 0)
  testthat::expect_lte(cm["Positive", "Negative"], 1)
})
