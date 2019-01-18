input1 <- list(mu_control = 10, mu_treatment = 8,
               sd_control = 0.8, sd_treatment = 1.2,
               N_total = 300, EndofStudy = 50,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270),
               N_impute = 20)

input2 <- list(mu_treatment = 8, sd_treatment = 1,
               N_total = 300, EndofStudy = 50,
               lambda = c(0.3, 1), lambda_time = c(25),
               interim_look = c(220, 270), h0 = 2,
               N_impute = 20)

context("")
test_that("The normal bayesCT is ", {
  set.seed(200)
  expect_equal(do.call(normalBACT, input1)$mu_treatment, 8)
  expect_equal(do.call(normalBACT, input2)$mu_treatment, 8)
  input1$sd_control <- -1.2
  expect_error(do.call(normalBACT, input1))

})
