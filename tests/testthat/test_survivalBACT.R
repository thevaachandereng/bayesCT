input1 <- list(hazard_treatment = 0.01, hazard_control = 0.02,
               N_total = 300,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270), EndofStudy = 50)


input2 <- list(hazard_treatment = c(0.01, 0.02), cutpoint = 25,
               EndofStudy = 50,N_impute = 10, N_total = 300,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270), alternative = "less")


context("")
test_that("The binomial bayesCT is ", {
  suppressWarnings(RNGversion("3.5.0"))
  expect_equal(do.call(survivalBACT, input1)$hazard_treatment, 0.01)
  expect_equal(do.call(survivalBACT, input2)$margin, 0.5)
  expect_equal(do.call(survivalBACT, input2)$prob_of_accepting_alternative, 0.95)
  input1$alternative <- "lessthan"
  expect_error(do.call(binomialBACT, input1))
})






