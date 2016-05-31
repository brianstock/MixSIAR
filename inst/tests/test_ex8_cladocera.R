# Test cladocera example
context("Ex script 8/9 (cladocera)")

test_that("Cladocera ex works",{
  mix.filename <- system.file("extdata", "cladocera_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("c14.0","c16.0","c16.1w9","c16.1w7","c16.2w4",
                                   "c16.3w3","c16.4w3","c17.0","c18.0","c18.1w9",
                                   "c18.1w7","c18.2w6","c18.3w6","c18.3w3","c18.4w3",
                                   "c18.5w3","c20.0","c22.0","c20.4w6","c20.5w3",
                                   "c22.6w3","BrFA"),
                       factors="id",
                       fac_random=FALSE,
                       fac_nested=FALSE,
                       cont_effects=NULL)
  source.filename <- system.file("extdata", "cladocera_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="means",
                             mix)
  discr.filename <- system.file("extdata", "cladocera_discrimination.csv", package = "MixSIAR")
  discr <- load_discr_data(filename=discr.filename, mix)

  model_filename <- "MixSIAR_model.txt"
  resid_err <- FALSE
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
