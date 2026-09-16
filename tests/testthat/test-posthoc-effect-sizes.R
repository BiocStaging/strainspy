testthat::test_that("comp_ani_diff_and_posthoc_test annotates top hits", {
  testthat::skip_if_not_installed("emmeans")

  d <- ss_test_data()
  fit <- ss_fit_zib()
  se_sub <- d$se_filtered[1:20, ]

  th <- top_hits(fit, coef = 2, alpha = 1)
  testthat::skip_if(nrow(th) == 0, "no top hits to post-hoc test")

  th <- head(th, 5)

  out <- suppressWarnings(
    comp_ani_diff_and_posthoc_test(
      se_sub, fit, th,
      nthreads = 1, progress = FALSE
    )
  )

  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_equal(nrow(out), nrow(th))
  testthat::expect_true("ANI_Difference" %in% colnames(out))
  testthat::expect_true("Comment" %in% colnames(out))
})
