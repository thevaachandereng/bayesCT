context("")
test_that("The binomial_outcome are", {
  expect_equal(binomial_outcome(p_control = 0.91, p_treatment = 0.83)$p_control, 0.91)
  expect_equal(binomial_outcome(p_control = 0.35, p_treatment = 0.32)$p_treatment, 0.32)
})

context("")
test_that("The beta non-informative priors are", {
  expect_equal(beta_prior(a0 = 0.5, b0 = 0.5)$prior, c(0.5, 0.5))
  expect_equal(beta_prior(a0 = 2, b0 = 2)$prior, c(2, 2))
})


context("")
test_that("The data binomial are", {
  expect_equal(data_binomial(treatment = binomialdata$treatment,
                             outcome   = binomialdata$outcome,
                             complete  = binomialdata$complete)$outcome,
               binomialdata$outcome)
})
