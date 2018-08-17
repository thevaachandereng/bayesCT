context("")
test_that("The binomial_outcome are", {
  expect_equal(binomial_outcome(p_control_true = 0.91, p_treatment_true = 0.83)$p_control, 0.93)
  expect_equal(binomial_outcome(p_control_true = 0.35, p_treatment_true = 0.32)$p_treatment, 0.93)
})
