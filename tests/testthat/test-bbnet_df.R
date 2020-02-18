

test_that("extract_stap_formula works", {
  expect_equal(sap_covs(extract_stap_data(formula = ~ sap(HFS))),
               c("HFS"))
  expect_identical(extract_stap_data(formula = ~(sap(HFS)|id))$group_indicator,c(1))
  expect_identical(extract_stap_data(formula = ~(sap(HFS)|id))$group_term,"id")
})

# test_that("extract_stap_lmer_formula works",{
#   
# })


# test_that("bbnet_df works", {
# 
# })
