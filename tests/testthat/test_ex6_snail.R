# Test snail example
context("Ex script 6/9 (snail)")

test_that("Snail ex works",{
  mix.filename <- system.file("extdata", "snail_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects=NULL)
  source.filename <- system.file("extdata", "snail_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="raw",
                             mix)
  discr.filename <- system.file("extdata", "snail_discrimination.csv", package = "MixSIAR")
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
