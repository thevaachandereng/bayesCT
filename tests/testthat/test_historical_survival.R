context("")
test_that("The historical_normal are", {
  expect_equal(historical_survival(time       = rexp(10, 1),
                                   treatment  = rep(1, 10),
                                   event      = rbinom(10, 1, 0.5))$treatment0, rep(1, 10))
  expect_equal(historical_survival(time               = rexp(10, 1),
                                   treatment          = rep(1, 10),
                                   event              = rbinom(10, 1, 0.5),
                                   discount_function  = "weibull")$discount_function, "weibull")
})
