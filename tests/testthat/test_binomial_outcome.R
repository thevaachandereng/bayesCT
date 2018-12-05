context("")
test_that("The binomial_outcome are", {
  expect_equal(binomial_outcome(p_control = 0.91, p_treatment = 0.83)$p_control, 0.91)
  expect_equal(binomial_outcome(p_control = 0.35, p_treatment = 0.32)$p_treatment, 0.32)
})
