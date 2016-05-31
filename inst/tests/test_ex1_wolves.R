# Test wolves example
context("Ex script 1/9 (wolves)")

test_that("Wolves ex works",{
  mix.filename <- system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=c("Region","Pack"),
                       fac_random=c(TRUE,TRUE),
                       fac_nested=c(FALSE,TRUE),
                       cont_effects=NULL)

  source.filename <- system.file("extdata", "wolves_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename, source_factors="Region",
                             conc_dep=FALSE, data_type="means", mix)

  discr.filename <- system.file("extdata", "wolves_discrimination.csv", package = "MixSIAR")
  discr <- load_discr_data(filename=discr.filename, mix)

  # plot_data(filename="isospace_plot",
  #           plot_save_pdf=FALSE,
  #           plot_save_png=FALSE,
  #           mix,source,discr)
  # calc_area(source=source,mix=mix,discr=discr)
  # plot_prior(alpha.prior=1,source)

  model_filename <- "MixSIAR_model.txt"
  resid_err <- TRUE
  process_err <- TRUE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

  run <- list(chainLength=3, burn=1, thin=1, chains=3, calcDIC=TRUE)
  invisible(capture.output(
    jags.1 <- run_model(run, mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)
  ))

  expect_is(jags.1,"rjags")
  file.remove(model_filename)
})
