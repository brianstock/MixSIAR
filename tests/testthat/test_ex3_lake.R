# Test lake example
context("Ex script 3/9 (lake)")

test_that("Lake ex works",{
  mix.filename <- system.file("extdata", "lake_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects="Secchi.Mixed")
  source.filename <- system.file("extdata", "lake_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="raw",
                             mix)
  discr.filename <- system.file("extdata", "lake_discrimination.csv", package = "MixSIAR")
  discr <- load_discr_data(filename=discr.filename, mix)

  model_filename <- "MixSIAR_model.txt"
  resid_err <- TRUE
  process_err <- FALSE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

  run <- list(chainLength=3, burn=1, thin=1, chains=3, calcDIC=TRUE)
  invisible(capture.output(
    jags.1 <- run_model(run, mix, source, discr, model_filename)
  ))

  expect_is(jags.1,"rjags")
  file.remove(model_filename)
})
