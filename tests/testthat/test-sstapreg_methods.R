SW <- suppressWarnings
bdf <- rbenvo::example_benvo
capture_output(
  lm1 <- SW(sstap_lm(BMI ~ sex  +sap(FFR), bdf,iter=500,chains=1))
  )
capture_output(
  glm1 <- SW(sstap_glm(cbind(NumObese,Num_NotObese) ~ sex + Centered_Scaled_Income + sap(HFS),
                       benvo = binomial_benvo,
                       family = binomial(),iter=500,chains=1))
)
capture_output(
  glm2 <- SW(sstap_glm(NumObese ~ sex + Centered_Scaled_Income + sap(HFS) + sap(FFR),
                       benvo = poisson_benvo,
                       family = poisson(),iter=500,chains=1))
)
capture_output(
  lmer1 <- SW(sstap_lmer(BMI ~ sex + year + stap(HFS) + (year|ID), benvo = complex_longitudinal,iter=200,chains=1))
)

test_that("Methods work on example models", {
  expect_equal(1000,nobs(lm1))
  expect_equal(1000,length(fitted(lm1)))
  expect_equal(c(12,12),dim(vcov(lm1)))
  expect_equal(300,nobs(glm1))
  expect_equal(300,nobs(glm2))
  expect_equal(c(15,2),dim(posterior_interval(lm1)))
  expect_equal(1171,nobs(lmer1))
})

# 
test_that("Sample glm with multiple predictors works",{
  
  expect_equal(1,length(glm1$specification$term))
  expect_equal(2,length(glm2$specification$term))
})

test_that("sstap_glmer methods work as intended",{
  expect_equal(1,length(ranef(lmer1)))
  expect_error(ranef(glm1))
  expect_equivalent(600,ngrps(lmer1))
})
