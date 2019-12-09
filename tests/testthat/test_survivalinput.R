context("")
test_that("The survival_outcome are", {
  expect_equal(survival_outcome(hazard_treatment = 0.01, hazard_control = 0.02)$hazard_control, 0.02)
  expect_equal(survival_outcome(hazard_treatment = c(0.01, 0.005), cutpoint = 25)$cutpoint, 25)
})

context("")
test_that("The non-informative gamma priors are", {
  expect_equal(gamma_prior(a0 = 0.1, b0 = 0.1)$prior, c(0.1, 0.1))
  expect_equal(gamma_prior(a0 = .2, b0 = .2)$prior, c(.2, .2))
})


context("")
test_that("The data_survival are", {
  expect_equal(data_survival(treatment = survivaldata$treatment,
                             time      = survivaldata$time,
                             event     = survivaldata$event)$time,
               survivaldata$time)
})
