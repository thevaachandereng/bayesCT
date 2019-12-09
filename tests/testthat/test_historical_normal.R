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
