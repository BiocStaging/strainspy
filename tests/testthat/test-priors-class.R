testthat::test_that("compute_eb_priors returns a valid strainspy_priors object", {
  testthat::skip_if_not_installed("fastglm")
  testthat::skip_if_not_installed("betareg")

  pri <- ss_priors_eb()

  testthat::expect_s4_class(pri, "strainspy_priors")
  testthat::expect_true(methods::validObject(pri))
  testthat::expect_identical(pri@method, "empirical")
  testthat::expect_s3_class(getPriorsDf(pri), "data.frame")
})

testthat::test_that("strainspy_priors show and accessor methods run", {
  testthat::skip_if_not_installed("fastglm")
  testthat::skip_if_not_installed("betareg")

  pri <- ss_priors_eb()

  testthat::expect_output(show(pri))

  fx <- extract_fixef_priors(pri)
  testthat::expect_s3_class(fx, "data.frame")
  testthat::expect_true(nrow(fx) > 0)
})

testthat::test_that("plot_prior_bootstrap draws for empirical priors", {
  testthat::skip_if_not_installed("fastglm")
  testthat::skip_if_not_installed("betareg")
  testthat::skip_if_not_installed("patchwork")

  pri <- ss_priors_eb()
  prior_names <- rownames(pri@boot_fixef)
  testthat::skip_if(length(prior_names) == 0, "no bootstrap draws available")

  p <- plot_prior_bootstrap(pri, prior_names[1])
  testthat::expect_false(is.null(p))

  testthat::expect_error(
    plot_prior_bootstrap(pri, "not_a_coefficient"),
    "not found in bootstrap matrices"
  )
})

testthat::test_that("define_priors covers preset and manual methods", {
  d <- ss_test_data()
  design <- stats::as.formula("~ Case_status + Age_at_collection")

  weak <- suppressWarnings(define_priors(d$se_filtered[1:10, ], design,
                                         method = "preset_weak"))
  testthat::expect_s4_class(weak, "strainspy_priors")

  strong <- suppressWarnings(define_priors(d$se_filtered[1:10, ], design,
                                           method = "preset_strong"))
  testthat::expect_s4_class(strong, "strainspy_priors")

  testthat::expect_error(
    define_priors(d$se_filtered[1:10, ], design, method = "manual"),
    "priors_df"
  )
})
