# Test "calc_area" function using wolves data
context("Calculate polygon area")

test_that("calc_area returns numeric",{
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

  expect_is(calc_area(source,mix,discr),"numeric")
})

test_that("calc_area throws error for != 2 iso",{
  mix.filename <- system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=c("Region","Pack"),
                       fac_random=c(TRUE,TRUE),
                       fac_nested=c(FALSE,TRUE),
                       cont_effects=NULL)
  mix$n.iso <- 3
  expect_error(calc_area(source,mix,discr))

  mix$n.iso <- 1
  expect_error(calc_area(source,mix,discr))
})
