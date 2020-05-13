# Test isopod example
context("Ex script 9/9 (isopod)")

test_that("Isopod ex works",{
  mix.filename <- system.file("extdata", "isopod_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("c16.4w3","c18.2w6","c18.3w3","c18.4w3","c20.4w6","c20.5w3","c22.5w3","c22.6w3"),
                       factors="Site",
                       fac_random=TRUE,
                       fac_nested=FALSE,
                       cont_effects=NULL)
  source.filename <- system.file("extdata", "isopod_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="means",
                             mix)
  discr.filename <- system.file("extdata", "isopod_discrimination.csv", package = "MixSIAR")
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
