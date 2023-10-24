test_that("writeFam() + readFam() recreates ped", {
  x1 = nuclearPed()
  y1 = writeFam(x1, famfile = tempfile(), verbose = F) |> readFam(verbose = F)
  expect_equal(x1, y1)

  # With marker
  x2 = x1 |> addMarker(geno = c(NA, "2/2", "1/3"), name = "M")
  y2 = writeFam(x2, famfile = tempfile(), verbose = F) |> readFam(verbose = F)
  expect_equal(x2, y2)

  # With empty marker w/ mutation model
  x3 = x2 |> addMarker(name = "M2", alleles = 1:2) |>
    setMutmod(markers = "M2", model = "equal", rate = 0.1)
  y3 = writeFam(x3, famfile = tempfile(), verbose = F) |> readFam(verbose = F)

  # x3 and y3 are not identical because y3 contains extra 0's in mut model
  y3 = setMutmod(y3, marker = "M2", rate2 = NULL, range = NULL, update = TRUE)

  expect_equal(x3, y3)

  # With stepwise model at both markers
  x4 = setMutmod(x3, model = "stepwise", rate = 0.01, rate2 = 1e-6, range = 0.1)
  y4 = writeFam(x4, famfile = tempfile(), verbose = F) |> readFam(verbose = F)
  expect_equal(x4, y4)
})
