test_that("unit conversions use correct molecular weights", {
  expect_equal(ngml_to_nM(1, "NTX"), 1000 / 341.4, tolerance = 1e-10)
  expect_equal(round(ngml_to_nM(1, "NTX"), 3), 2.929)
  expect_equal(nM_to_ngml(ngml_to_nM(7.3, "BN"), "BN"), 7.3)
  expect_equal(mg_to_nmol(4.5, "NTX"), 4.5e6 / 341.4)
})

test_that("logit and inv_logit are inverses", {
  expect_equal(inv_logit(logit(c(0.01, 0.5, 0.97))), c(0.01, 0.5, 0.97))
})
