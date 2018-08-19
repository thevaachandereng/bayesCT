context("")
test_that("The hypothesis are", {
  expect_equal(hypothesis(0, 0.05, 0.95, 0.90)$h0, 0)
  expect_equal(hypothesis(0, 0.05, 0.95, 0.90)$futility_prob, 0.05)
  expect_equal(hypothesis(0, 0.05, 0.95, 0.90)$prob_ha, 0.95)
  expect_equal(hypothesis(0, 0.05, 0.95, 0.90)$expected_success_prob, 0.90)
})
