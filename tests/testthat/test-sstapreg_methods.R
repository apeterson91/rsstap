SW <- suppressWarnings

# capture_output(
#   lm1 <- SW(sstap_lm(BMI ~ sex  +sap(FFR), example_benvo,iter=500,chains=1)),
#   glm1 <- SW(sstap_glm(cbind(NumObese,Num_NotObese) ~ sex + Centered_Scaled_Income + sap(HFS),
#                     benvo = binomial_benvo,
#                     family = binomial(),iter=500,chains=1)),
#   glm2 <- SW(sstap_glm(NumObese ~ sex + Centered_Scaled_Income + sap(HFS) + sap(FFR),
#                        benvo = poisson_benvo,
#                        family = poisson(),iter=500,chains=1))
# )
# 
# test_that("Methods work on example models", {
#   expect_equal(1000,nobs(lm1))
#   expect_equal(1000,length(fitted(lm1)))
#   expect_equal(c(2,2),dim(vcov(lm1)))
# })
# 
# test_that("Sample glm with multiple predictors works",{
#   expect_equal(2,length(glm$stap_terms))
# })
