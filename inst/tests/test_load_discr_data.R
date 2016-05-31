# Test "load_discr_data" function using wolves data
context("Load discrimination data")

test_that("Discrimination data matches source data",{
  # Correctly loaded wolves data
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

  # Not test
  # colnames(DISCR) contains mix$MU_names and mix$SIG_names (e.g. 'Meand13C')

  # check that discr rownames match source_names
  # throws error message in 'plot_data.R'
  expect_identical(source$source_names,rownames(discr$mu))
  expect_identical(source$source_names,rownames(discr$sig2))
})
