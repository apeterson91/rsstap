bdf <- rbenvo::longitudinal_HFS
test_that("checks error as appropriate", {
  expect_error(sstap_lmer(BMI ~ sex + (1|id) ,
                          benvo=bdf),
               regexp = "No covariates designated as")
})
