# Test killer whale example
context("Ex script 5/9 (killer whale)")

test_that("Killer whale ex works",{
  mix.filename <- system.file("extdata", "killerwhale_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects=NULL)
  source.filename <- system.file("extdata", "killerwhale_sources.csv", package = "MixSIAR")
  source <- load_source_data(filename=source.filename,
                             source_factors=NULL,
                             conc_dep=FALSE,
                             data_type="means",
                             mix)
  discr.filename <- system.file("extdata", "killerwhale_discrimination.csv", package = "MixSIAR")
  discr <- load_discr_data(filename=discr.filename, mix)

  # Test uninformative prior
  model_filename <- "MixSIAR_model_kw_uninf.txt"   # Name of the JAGS model file
  resid_err <- TRUE
  process_err <- TRUE
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)

  run <- list(chainLength=5, burn=2, thin=1, chains=3, calcDIC=TRUE)
  invisible(capture.output(
    jags.uninf <- run_model(run,mix,source,discr,model_filename,alpha.prior=1)
  ))
  expect_is(jags.uninf,"rjags")
  file.remove(model_filename)

  # Test informative prior version
  kw.alpha <- c(10,1,0,0,3)
  kw.alpha <- kw.alpha*length(kw.alpha)/sum(kw.alpha)
  kw.alpha[which(kw.alpha==0)] <- 0.01
  model_filename <- "MixSIAR_model_kw_inf.txt"   # Name of the JAGS model file
  write_JAGS_model(model_filename, resid_err, process_err, mix, source)
  invisible(capture.output(
    jags.inf <- run_model(run=run,mix,source,discr,model_filename,alpha.prior=kw.alpha)
  ))
  expect_is(jags.inf,"rjags")
  file.remove(model_filename)
})
