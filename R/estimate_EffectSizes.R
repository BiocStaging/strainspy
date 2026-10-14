#' Compute the mean ANI differences for contigs with an identity effect and 
#' perform posthoc testing to establish support for each hit.
#'
#' In the zero-inflated beta model, beta effects represent changes in
#' strain identity and are modeled on the logit scale. While this scale
#' is effective for distinguishing phenotype groups statistically, the
#' resulting effect sizes are not directly interpretable as differences
#' in mean ANI. The same applies to effect sizes from the ordinal beta model.
#'
#' This function converts logit-scale beta effects into response-scale
#' mean ANI differences using marginal means and pairwise contrasts
#' (requires `emmeans`). It therefore provides an interpretable estimate of the
#' magnitude of strain-level differences between phenotype groups, along
#' with post hoc quality checks for poorly supported beta signals. 
#'
#' Posthoc testing uses two approaches. First, it counts the number of non-zero 
#' observations supporting each beta effect. Second, it compares the estimated 
#' mean ANI difference against a preset threshold. A comment is added to each 
#' hit depending on the outcome.
#' 
#' By default, this function uses the formula from `getDesign(fit)`. However, if the
#' formula contains multiple covariates, this can cause emmeans to behave 
#' unexpectedly. For example, if `getDesign(fit)` is `~ disease + age + sex + country + (1|study)`,
#' this can cause this function to throw errors or warnings.  Using a `simplified_formula = ~ disease`
#' to override can be a useful approach.
#'
#' This function is designed to directly work on the `tibble` returned from 
#' `top_hits()`.
#' 
#' @param se SummarizedExperiment. A `SummarizedExperiment` object containing 
#' assay data and metadata
#' @param fit A fitted `strainspy_fit`model
#' @param th The table returned by `top_hits()` (required).
#' @param simplified_formula A formula to use instead of `getDesign(fit)`, useful
#' when the original design was too complex and causes emmeans to fail during 
#' ANI estimation (default NULL).
#' @param beta_min_nz Minimum proportion of non-zero values per phenotype group for
#' a beta hit to be considered well supported. Hits where that proportion falls below
#' this threshold in at least one group are flagged `too few nonzero`. Default 0.1.
#' Hits supported by 3 or fewer observations are always flagged `<3 nonzero`.
#' @param beta_min_ani_diff Minimum absolute mean ANI difference required  for 
#' a beta hit to be considered meaningful. Effect sizes smaller than this are 
#' flagged: `small effect`. Default = 1.5e-2, i.e. 1.5 percentage points of ANI. 
#' @param nthreads Number of threads for parallel processing. Default `1`.
#' @param reorder_hits Logical. Return `th` after reordering hits to prioritise
#' the best-supported ones. Default `FALSE`.
#' @param BPPARAM Optional BiocParallelParam object.
#' @param progress Logical. If `TRUE` (default), progress bars and progress messages are printed. Set to `FALSE` to silence them.
#'
#' @return A `data.frame` similar to `th`, with added columns giving details on:
#' Component (Beta, ZI or both), Mean ANI_Difference, Contrast, Non-Zero Counts 
#' of contrast components, and a general `Comment` on each hit. The comment will 
#' be `NA` if the hit passes all posthoc tests.
#' @examples
#' if (requireNamespace("emmeans", quietly = TRUE)) {
#'   meta <- read.csv(system.file("extdata", "example_metadata.csv.gz", 
#'   package = "strainspy"))
#'   meta$Case_status <- factor(meta$Case_status)
#'   se <- read_sylph(system.file("extdata", "example_sylph_profile.tsv.gz", 
#'   package = "strainspy"), meta_data = meta)
#'   se <- filter_by_presence(se, min_nonzero = 30)[1:10, ]
#'   fit <- glmZiBFit(se, ~ Case_status + Age_at_collection, nthreads = 1)
#'   th <- top_hits(fit, coef = 2, alpha = 1)
#'   posthoc <- comp_ani_diff_and_posthoc_test(se, fit, th, nthreads = 1)
#'   head(posthoc$ANI_Difference)
#' }
#' 
#' @export
#' @importFrom BiocParallel bplapply SnowParam SerialParam
comp_ani_diff_and_posthoc_test <- function(se, fit, th = NULL,
                                          simplified_formula = NULL,
                                          beta_min_nz = 0.1, 
                                          beta_min_ani_diff = 1.5e-2, 
                                          nthreads = 1L, reorder_hits = FALSE,
                                          BPPARAM=NULL, progress = TRUE){
  
  check_progress(progress)
  
  require_suggested("emmeans", "post-hoc effect size estimation")
  
  # check thresholds
  if (!is.numeric(beta_min_nz) || beta_min_nz <= 0 || beta_min_nz >= 1)
    stop("`beta_min_nz` must be between 0 and 1.")
  
  if (!is.numeric(beta_min_ani_diff) || beta_min_ani_diff < 0)
    stop("`beta_min_ani_diff` must be >= 0.")
  
  if (!inherits(se, "SummarizedExperiment"))
    stop("`se` must be a SummarizedExperiment.")
  
  if (!inherits(fit, "strainspy_fit"))
    stop("`fit` must be a strainspy_fit object.")
  
  if(is.null(getVcov(fit)))
    stop("`fit` must contain variance-covariance data.")
  
  
  
  # Annotate the top hits
  th <- annotate_effect(th = th, pval_thresh = attr(th, 'alpha'))
  
  if(!"beta" %in% th$component){ # no beta hits
    if(!"both" %in% th$component) {
      stop("No identity difference hits found in the provided th: ",
           "running `comp_ani_diff_and_posthoc_test()` is not required.")
    }
  }
  
  top_hit_contigs <- th$Contig_name[th$component !="zi"]
  
  # Could be a wrong se are the contigs present?
  missing_contigs <- setdiff(top_hit_contigs, rownames(se))
  if (length(missing_contigs) > 0) {
    missing_txt <- paste(missing_contigs, collapse = ", ")
    stop("The following contigs are not present in `se`: ", missing_txt)
  } else {
    coef_idx <- match(top_hit_contigs, getContigNames(fit))
  }
  
  # pick col_data from se
  col_data <- SummarizedExperiment::colData(se)
  
  # stick with scaling done during the original fit
  if (isScaleContinuous(fit)){ 
    col_data <- scale_numeric_columns(col_data)
  }
  
  
  # Set up parallel infrastructure
  BPPARAM <- setup_bpparam(nthreads, BPPARAM, progress, tasks = length(top_hit_contigs))
  
  
  op_list <- BiocParallel::bplapply(
    coef_idx,
    process_contig,
    fit = fit,
    form = simplified_formula,
    se = se,
    col_data = col_data,
    pheno = attr(th, 'phenotype_coef'),
    BPPARAM = BPPARAM
  )
  
  op <- dplyr::bind_rows(op_list)
  
  op <- dplyr::left_join(th, op, by = "Contig_name")
  
  op <- posthoc_test(op, beta_min_nz, beta_min_ani_diff)
  
  if(reorder_hits){
    # Sort the hits in a more helpful way? Mixing beta + Zi is a bit of a headache
    op <- op %>%
      dplyr::arrange(
        dplyr::desc(is.na(fail_reasons)),           # valid before flagged (beta only)
        component,
        component_p
      )
  }
  
  
  return(data.frame(Contig_name = op$Contig_name,
                    Genome_file = op$Genome_file,
                    coefficient = op$coefficient,
                    std_error = op$std_error,
                    p_value = op$p_value,
                    p_adjust = op$p_adjust,
                    zi_coefficient = op$zi_coefficient,
                    zi_std_error = op$zi_std_error,
                    zi_p_value = op$zi_p_value,
                    zi_p_adjust = op$zi_p_adjust,
                    hit_component = op$component,
                    ANI_Difference = op$estimate,
                    contrast = op$contrast,
                    ref_nz_obs = op$ref_nz,
                    ref_tot_obs = op$ref_total,
                    comp_nz_obs = op$comp_nz,
                    comp_tot_obs = op$comp_total,
                    Comment = op$fail_reasons
  ))
  
}


process_contig <- function(
    i, # row_index
    fit,
    form = NULL,
    se,
    col_data,
    pheno
) {
  
  beta1  <- unlist(getCoefficients(fit)[i,])
  # Intercept name is sometimes off (?bug)
  names(beta1) <- sub(".*[Ii]ntercept.*", "(Intercept)", names(beta1))
  lvls <- c("(Intercept)", names(beta1)[pheno])
  beta1 <- beta1[lvls]
  
  if(is.null(form)){
    combined_formula <- as.formula(paste(c("Value", as.character(getDesign(fit))),
                                                        collapse = " "))
  } else{
    combined_formula <- as.formula(paste(c("Value", as.character(form)),
                                        collapse = " "))
  }
  
  
  
  formula1 <- nobars_(combined_formula)
  
  V1 <- getVcov(fit)[[i]][lvls, lvls]
  
  cd <- col_data
  cd$Value <- offset_ANI(
    as.vector(SummarizedExperiment::assay(se)[i, ]) / 100
  )
  cd_df <- as.data.frame(cd)
  
  # subset columns to only contain variables in formula
  mf1    <- cd_df[, all.vars(combined_formula)]
  # subset to only keep relavant rows
  mdlmx <- model.matrix(formula1, cd_df)
  
  ref_rows <- which(mdlmx[, "(Intercept)"] == 1 & mdlmx[,lvls[2]] == 0)
  comp_rows <- which(mdlmx[, lvls[2]] == 1)
  
  # find out support
  s1 <- quickly_find_support(ref_rows, cd_df)
  s2 <- quickly_find_support(comp_rows, cd_df)
  
  picked_rows <- sort(c(ref_rows, comp_rows))
  
  mf1 <- mf1[picked_rows, ]
  
  mf1 <- remove_bad_rows(mf1) # qdrg doesn't like bad rows, remove them and warn!
  
  rg <- emmeans::qdrg(
    formula = formula1,
    data = mf1,
    coef = beta1,
    vcov = V1,
    df = Inf,
    family = getFamily(fit)
  )
  
  emm1 <- emmeans::emmeans(rg, stats::delete.response(stats::terms(formula1)))
  res <- data.frame(emmeans::contrast(emm1, method = "trt.vs.ctrl"))
  pck_row <- which(!is.nan(res$p.value))
  
  if(length(pck_row) > 1){
    if(length(unique(res$p.value[pck_row]))){
      pck_row <- pck_row[1] # If there are multiple rows (happens when several covariates are present)
    }
  }
  
  if(length(pck_row) == 1){
    data.frame(Contig_name = rownames(se)[i], ref_nz = s1[1], ref_total = s1[2], comp_nz = s2[1], comp_total = s2[2], res[pck_row,])
  } else {
    data.frame(Contig_name = rownames(se)[i], ref_nz = s1[1], ref_total = s1[2], comp_nz = s2[1], comp_total = s2[2], data.frame(
      contrast = NA, estimate = NA, SE = NA, df = NA, z.ratio = NA, p.value = NA))
  }
  
}


quickly_find_support <- function(rows, data) {
  vals <- data$Value[rows]
  c(
    nonzero = sum(vals != 0),
    total   = length(vals)
  )
}

remove_bad_rows <- function(df) {
  if (!is.data.frame(df)) stop("Input must be a data.frame")
  
  # Identify rows with NA, NULL, or Inf
  bad_rows <- apply(df, 1, function(row) {
    any(is.na(row) | is.infinite(row) | vapply(row, is.null, logical(1)))
  })
  
  n_removed <- sum(bad_rows)
  
  if (n_removed > 0) {
    warning(sprintf("Removed %d rows containing NA, NULL, or Inf", n_removed))
  }
  
  # Return cleaned data.frame
  df[!bad_rows, , drop = FALSE]
}

annotate_effect <- function(top_hit_contigs = NULL, 
                            phenotype_coef = NULL, 
                            ZB_fit = NULL, 
                            pval_thresh = 0.05, 
                            th = NULL) {
  
  if (is.null(th)) {
    th <- top_hits(fit = ZB_fit, 
                   coef = phenotype_coef, 
                   alpha = 1, 
                   method = "BH") %>%
      dplyr::filter(Contig_name %in% top_hit_contigs)
  }
  
  th <- th %>%
    dplyr::mutate(
      component = dplyr::case_when(
        p_adjust < pval_thresh & zi_p_adjust < pval_thresh ~ "both",
        p_adjust < pval_thresh ~ "beta",
        zi_p_adjust < pval_thresh ~ "zi",
        TRUE ~ "none"
      ),
      component_p = dplyr::case_when(
        component == "beta" ~ p_adjust,
        component == "zi" ~ zi_p_adjust,
        component == "both" ~ pmin(p_adjust, zi_p_adjust),
        TRUE ~ NA_real_
      )
    )
  
  return(th)
}

posthoc_test <- function(op, beta_min_nz_ratio, beta_min_ani_diff){
  op <- op %>% dplyr::rowwise() %>%
    dplyr::mutate(
      fail_reasons = paste(
        c(
          if (component %in% c("beta", "both") & min(c(ref_nz, comp_nz)) <= 3) "<3 nonzero" else NULL,
          if (component %in% c("beta", "both") & min(c(ref_nz/ref_total , comp_nz/comp_total)) <= beta_min_nz_ratio) "too few nonzero" else NULL,
          if (component %in% c("beta", "both") & abs(estimate) < beta_min_ani_diff) "small effect" else NULL
        ),
        collapse = "; "
      ),
      fail_reasons = ifelse(fail_reasons == "", NA_character_, fail_reasons)
    ) %>%
    dplyr::ungroup()
  
  return(op)
}
