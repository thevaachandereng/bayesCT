input1 <- list(p_control = 0.12, p_treatment = 0.08,
               y0_treatment = 8, N0_treatment = 90, N_impute = 10,
               y0_control = 13, N0_control = 95, N_total = 300,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270), EndofStudy = 50)


input2 <- list(p_treatment = 0.09, EndofStudy = 50,
               y0_treatment = 8, N0_treatment = 90,
               N_impute = 10, N_total = 300,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270), alternative = "less")

input3 <- list(p_control = 0.12, p_treatment = 0.08, N_impute = 10,
               N_total = 500, lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(420, 470), EndofStudy = 50, alternative = "two-sided")

context("")
test_that("The binomial bayesCT is ", {
  suppressWarnings(RNGversion("3.5.0"))
  expect_equal(do.call(binomialBACT, input1)$p_treatment, 0.08)
  expect_equal(do.call(binomialBACT, input1)$p_control, 0.12)
  expect_equal(do.call(binomialBACT, input2)$p_treatment, 0.09)
  expect_equal(do.call(binomialBACT, input3)$prob_of_accepting_alternative, 0.95)
  input1$p_control <- 1.2
  expect_error(do.call(binomialBACT, input1))
  input1$p_control <- 0.5
  input1$alternative <- "lessthan"
  expect_error(do.call(binomialBACT, input1))
})






