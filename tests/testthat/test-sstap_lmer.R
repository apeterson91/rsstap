test_that("checks error as appropriate", {
  expect_error(sstap_lmer(BMI ~ sex + (1|ID) ,
                          benvo=rbenvo::longitudinal_HFS),
               regexp = "No covariates designated as")
})