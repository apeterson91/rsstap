bdf <- rbenvo::longitudinal_HFS
test_that("sstapspec object and methods work correctly", {
  expect_equal(BMI ~ sex,get_sstapspec(BMI ~ sex + stap(HFS),
                                       benvo = bdf)$stapless_formula)
  expect_equivalent(0L,has_bw(get_sstapspec(BMI ~ sex + stap(HFS),
                                                benvo = bdf)))
  expect_equivalent(1,has_bw(get_sstapspec(BMI ~ sex + stap_bw(HFS),
                                           benvo = bdf)))
  expect_error(get_sstapspec(BMI ~ sex,
                             benvo=bdf),regexp="No covariates")
})
