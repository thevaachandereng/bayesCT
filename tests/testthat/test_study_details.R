context("")
test_that("The study details are", {
  expect_equal(study_details(total_sample_size = 100, study_period = 50,
                             interim_look = c(80, 90))$N_total, 100)
  expect_equal(study_details(total_sample_size = 100, study_period = 50,
                             interim_look = c(80, 90))$interim_look, c(80, 90))
  expect_equal(study_details(total_sample_size = 100, study_period = 50,
                             interim_look = c(80, 90))$EndofStudy, 50)
})
