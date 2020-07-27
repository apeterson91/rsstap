test_that("Extract Stap Formula works", {
  expect_equal(BMI ~ sex,get_stapless_formula(BMI ~ sex + stap(HFS))$stapless_formula)
  expect_equivalent(c("0"),get_stapless_formula(BMI ~ sex + stap(HFS))$stap_mat[,3])
  expect_equivalent(c("1"),get_stapless_formula(BMI ~ sex + stap_bw(HFS))$stap_mat[,3])
  expect_error(get_stapless_formula(BMI ~ sex)$stapless_formula,regexp="No covariates")
})
