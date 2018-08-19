context("")
test_that("The historical_binomial are", {
  expect_equal(historical_binomial(y0_treatment = 10, y0_control = 12,
                                   N0_treatment = 29, N0_control = 25)$N0_control, 25)
  expect_equal(historical_binomial(y0_treatment = 10, y0_control = 12,
                                   N0_treatment = 29, N0_control = 25)$discount_function, "identity")
})
