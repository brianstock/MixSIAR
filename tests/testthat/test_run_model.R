# Test "run_model" function using wolves data
context("Run JAGS model")

test_that("Error messages work",{
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

  model_filename <- "MixSIAR_model.txt"
  resid_err <- TRUE
  process_err <- TRUE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

  run <- list(chainLength=3, burn=1, thin=1, chains=3, calcDIC=TRUE)

  # process and resid error cannot both be false
  expect_error(run_model(run, mix, source, discr, model_filename,
                        alpha.prior = 1, resid_err=FALSE, process_err=FALSE))

  # if mix$N==1, must choose process error only (MixSIR)
  mix$N <- 1
  expect_error(run_model(run, mix, source, discr, model_filename,
                         alpha.prior = 1, resid_err=TRUE, process_err=FALSE))
  expect_error(run_model(run, mix, source, discr, model_filename,
                         alpha.prior = 1, resid_err=TRUE, process_err=TRUE))

  # Error checks on alpha prior
  mix$N <- 66 # put N back to what it should be
  expect_error(run_model(run, mix, source, discr, model_filename,
                         alpha.prior = "blah", resid_err=TRUE, process_err=TRUE))
  expect_error(run_model(run, mix, source, discr, model_filename,
                         alpha.prior = 1:5, resid_err=TRUE, process_err=TRUE))
  expect_error(run_model(run, mix, source, discr, model_filename,
                         alpha.prior = 0:2, resid_err=TRUE, process_err=TRUE))

  # # cannot set informative prior on model with fixed effect
  # mix$n.fe <- 1
  # expect_error(run_model(run, mix, source, discr, model_filename,
  #                        alpha.prior = 1:3, resid_err=TRUE, process_err=TRUE))
  file.remove(model_filename)
})
