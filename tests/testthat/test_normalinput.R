context("")
test_that("The binomial_outcome are", {
  expect_equal(normal_outcome(mu_control = 0.91, mu_treatment = 0.88,
                              sd_control = 0.05, sd_treatment = 0.07)$sd_control, 0.05)
  expect_equal(normal_outcome(mu_control = 3.8, mu_treatment = 3.3,
                              sd_control = 1.2, sd_treatment = 1.6)$mu_control, 3.8)
})

context("")
test_that("The historical_normal are", {
  expect_equal(historical_normal(mu0_treatment = 15, sd0_treatment = 2, N0_treatment = 10,
                                 mu0_control = 17, sd0_control = 3, N0_control = 20,
                                 discount_function = "weibull", alpha_max = 1,
                                 fix_alpha = "FALSE", weibull_scale = 0.135, weibull_shape = 3,
                                 method = "mc")$sd0_control, 3)
  expect_equal(historical_normal(mu0_treatment = 15, sd0_treatment = 2, N0_treatment = 10,
                                 mu0_control = 17, sd0_control = 3, N0_control = 20)$discount_function, "identity")
})


context("")
test_that("The data binomial are", {
  expect_equal(data_normal(treatment = normaldata$treatment,
                           outcome   = normaldata$outcome,
                           complete  = normaldata$complete)$outcome,
               normaldata$outcome)
})
