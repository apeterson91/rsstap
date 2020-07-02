test_that("Extract Stap Formula works", {
  expect_equal(BMI ~ sex,get_stapless_formula(BMI ~ sex + stap(HFS))$stapless_formula)
  expect_error(get_stapless_formula(BMI ~ sex)$stapless_formula,regexp="No covariates")
})
