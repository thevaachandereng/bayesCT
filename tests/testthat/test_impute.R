context("")
test_that("The impute are", {
  expect_equal(impute(no_of_impute = 100, number_mcmc = 10000)$N_impute, 100)
  expect_equal(impute(no_of_impute = 100, number_mcmc = 1000)$number_mcmc, 1000)
})

