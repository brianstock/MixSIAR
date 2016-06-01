# Test "load_source_data" function using wolves data
context("Load source data")

test_that("Error messages work",{
  # Correctly loaded wolves mixture data
  mix.filename <- system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")
  mix <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=c("Region","Pack"),
                       fac_random=c(TRUE,TRUE),
                       fac_nested=c(FALSE,TRUE),
                       cont_effects=NULL)

  source.filename <- system.file("extdata", "wolves_sources.csv", package = "MixSIAR")

  # 1. cannot have > 1 source factor
  expect_error(load_source_data(filename=source.filename,
                                source_factors=c("Region","Pack"),
                                conc_dep=FALSE,
                                data_type="means", mix))

  # 2. source factor must be in colnames(SOURCE)
  expect_error(load_source_data(filename=source.filename,
                                source_factors=c("blah"),
                                conc_dep=FALSE,
                                data_type="means", mix))

  # 3. source factor must be in mix$factors
  # load mix data again with NO FACTORS
  mix2 <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=NULL,
                       fac_random=NULL,
                       fac_nested=NULL,
                       cont_effects=NULL)
  expect_error(load_source_data(filename=source.filename,
                                source_factors=c("Region"),
                                conc_dep=FALSE,
                                data_type="means", mix2))

  # 4. If conc_dep, must be of form 'Concd13C'
  expect_error(load_source_data(filename=source.filename,
                                source_factors=c("Region"),
                                conc_dep=TRUE,
                                data_type="means", mix))

  # 5. If raw source data, mix$iso_names must be in colnames(SOURCE)
  #    not tested here - would need to create another .csv data file
  #    lines 121-124 of 'load_source_data.R'

  # 6. If means + SD source data, 'Mean' and 'SD' + mix$iso_names must be
  #    in colnames(SOURCE). Not tested here as above.

  # 7. If means + SD source data, 'n' must be in colnames(SOURCE).
  #    Not tested here as above.

  # 8. If source SD = 0, throw error. Not tested here as above.
})
