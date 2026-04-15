test_that("hte3_example_data returns the expected columns", {
  data <- hte3_example_data(n = 25, seed = 1)

  expect_s3_class(data, "data.table")
  expect_true(all(c("W1", "W2", "W3", "A", "Y", "Y_binary", "tau", "pi1", "mu0", "mu1") %in% names(data)))
  expect_equal(nrow(data), 25)
})

test_that("truncate_propensity keeps scores in bounds", {
  scores <- c(0, 0.01, 0.5, 0.99, 1)
  treatment <- c(0, 0, 1, 1, 0)

  bounded <- truncate_propensity(scores, treatment, treatment_level = 1, truncation_method = "deterministic")

  expect_true(all(bounded > 0))
  expect_true(all(bounded < 1))
})
