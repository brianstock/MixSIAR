# Test "load_mix_data" function using wolves data
context("Load mixture data")

test_that("Error messages work",{
	mix.filename <- system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")

	# Number of isotopes must be > 0
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=NULL,
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# length of 'fac_random' must match length of 'factors'
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# Cannot have > 2 'factors'
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack","blah"),
                     			fac_random=c(TRUE,TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# length of fac_nested must match length of factors
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# if 2 factors, length of fac_nested must = 2
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE),
                     			cont_effects=NULL))
	# if 2 factors, both cannot be nested within each other
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(TRUE,TRUE),
                     			cont_effects=NULL))
	# cannot have > 1 continuous effect
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=c("blah1","blah2")))
	# check that iso_names must be in colnames(mix.filename)
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("D13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# check that factors must be in colnames(mix.filename)
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects=NULL))
	# check that cont_effects must be in colnames(mix.filename)
	expect_error(load_mix_data(filename=mix.filename,
                     			iso_names=c("d13C","d15N"),
                     			factors=c("Region","Pack"),
                     			fac_random=c(TRUE,TRUE),
                     			fac_nested=c(FALSE,TRUE),
                     			cont_effects="blah"))
})
