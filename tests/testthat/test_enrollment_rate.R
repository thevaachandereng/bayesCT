context("")
test_that("The enrollment rate are", {
  expect_equal(enrollment_rate(lambda = c(0.3, 0.9), time = 10)$lambda, c(0.3, 0.9))
  expect_equal(enrollment_rate(lambda = c(0.3, 0.9), time = 10)$lambda_time, 10)
})
