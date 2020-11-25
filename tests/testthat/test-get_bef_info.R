test_that("BEF Info", {
  expect_equivalent(c("sap","FFR",")","0"),
               get_bef_info(c("sap","tap","stap","gap"),
                            BMI ~ sex + sap(FFR)))
  expect_equivalent(c("tap","FFR",")","0"),
                    get_bef_info(c("sap","tap","stap","gap"),
                                 BMI ~ sex + tap(FFR)))
  expect_equivalent(c("tap","FFR",")","1"),
                    get_bef_info(c("sap","tap","stap","gap"),
                                 BMI ~ sex + tap_bw(FFR)))
  expect_equal(NULL,get_bef_info(c("sap","tap","stap","gap"),
                                 BMI ~ sex ))
})
