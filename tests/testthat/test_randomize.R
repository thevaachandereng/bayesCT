context("")
test_that("The randomize are", {
  expect_equal(randomize(block_size = 10, randomization_ratio = c(1, 1))$block, 10)
  expect_equal(randomize(block_size = 10, randomization_ratio = c(1, 1))$rand_ratio, c(1, 1))
})
