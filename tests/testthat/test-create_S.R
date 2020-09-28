ref_one <- mgcv::jagam(formula = id ~ -1 + s(Distance,bs='ps'), family = gaussian(),
                       data = rbenvo::joinvo(rbenvo::FFbenvo,"FFR","Distance",
                                             NA_to_zero = TRUE),
                       file = tempfile(fileext = ".jags"),
                       offset = NULL,
                       centred = FALSE,
                       diagonalize = FALSE)

ref_two <- mgcv::jagam(formula = id ~ -1 + s(Time,bs='ps'), family = gaussian(),
                       data = rbenvo::joinvo(rbenvo::longitudinal_HFS,"HFS","Time",
                                             NA_to_zero = TRUE),
                       file = tempfile(fileext = ".jags"),
                       offset = NULL,
                       centred = FALSE,
                       diagonalize = FALSE)

ref_three <- mgcv::jagam(formula = id ~ -1 + t2(Distance,Time,bs='ps'), family = gaussian(),
                       data = rbenvo::joinvo(rbenvo::longitudinal_HFS,"HFS","Distance-Time",
                                             NA_to_zero = TRUE),
                       file = tempfile(fileext = ".jags"),
                       offset = NULL,
                       centred = FALSE,
                       diagonalize = FALSE)

test_that("create_S works", {
  ## No interaction or decompositions
  expect_equal(c(10,20),dim(create_S(ref_one,0)))
  expect_equal(c(10,20),dim(create_S(ref_two,0)))
  expect_equal(c(25,100),dim(create_S(ref_three,0)))
  ## between_within decompositions
  expect_equal(list(c(10,20),c(10,20)),lapply(create_S(ref_one,1),dim))
  expect_equal(list(c(10,20),c(10,20)),lapply(create_S(ref_two,1),dim))
  expect_equal(list(c(25,100),c(25,100)),lapply(create_S(ref_three,1),dim))
  ## Interactions

})
