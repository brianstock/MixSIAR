# Test geese example
context("Ex script 2/9 (geese)")

test_that("Geese ex works",{
  mix.filename <- system.file("extdata", "geese_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors="Group",
                       fac_random=FALSE,
                       fac_nested=FALSE,
                       cont_effects=NULL)
  source.filename <- system.file("extdata", "geese_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=TRUE,
                             data_type="means",
                             mix)
  discr.filename <- system.file("extdata", "geese_discrimination.csv", package = "MixSIAR")
  discr <- load_discr_data(filename=discr.filename, mix)

  model_filename <- "MixSIAR_model.txt"
  resid_err <- TRUE
  process_err <- FALSE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

  run <- list(chainLength=3, burn=1, thin=1, chains=3, calcDIC=TRUE)
  invisible(capture.output(
    jags.1 <- run_model(run, mix, source, discr, model_filename,
                      alpha.prior = 1, resid_err, process_err)
  ))

  expect_is(jags.1,"rjags")
  file.remove(model_filename)
})
