input <- list(p_control = 0.12, p_treatment = 0.10,
              y0_treatment = 8, N0_treatment = 90, N_impute = 10,
              y0_control = 13, N0_control = 95, N_total = 300,
              lambda = c(0.3, 1), lambda_time = c(25),
              interim_look = c(220, 270), EndofStudy = 50)


context("")
test_that("The binomial bayesCT is ", {
  set.seed(200)
  expect_equal(do.call(binomialBACT, input)$p_treatment, 0.10)
  input$p_control <- 1.2
  expect_error(do.call(binomialBACT, input))
})






