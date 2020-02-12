library(testthat)
library(DarkDiv)

test_that("Beals gives same results as vegan", {
  require(vegan)
  data(dune)
  BealsDD <- DarkDiv(x = dune, method = "RawBeals")
  BealsVegan <- vegan::beals(x = dune, type = 0, include=F)
  expect_equal(BealsDD$AllProbs, BealsVegan)
})

test_that("Using indication matrix directly yields same results as reference", {
  require(vegan)
  data(dune)
  #RawBeals
  ddBeals <- DarkDiv(x = dune, method = "RawBeals")
  ddBeals2 <- DarkDiv(x = dune, r = ddBeals$indication,  method = "RawBeals")
  expect_equal(ddBeals, ddBeals2)
  #Favorability
  ddFav <- DarkDiv(x = dune, method = "Favorability")
  ddFav2 <- DarkDiv(x = dune, r = ddFav$indication,  method = "Favorability")
  expect_equal(ddFav, ddFav2)
  #Hypergeometric
  ddHyp <- DarkDiv(x = dune, method = "Hypergeometric")
  ddHyp2 <- DarkDiv(x = dune, r = ddHyp$indication,  method = "Hypergeometric")
  expect_equal(ddHyp, ddHyp2)
})

test_that("When removeAbsent = F, names of species and sites are unchanged, even if some species are removed", {
  require(vegan)
  data(dune)
  dune2 <- dune
  dune2[ , 2] <- rep(0, nrow(dune2))
  #RawBeals
  ddHypFull <- DarkDiv(x = dune, method = "Hypergeometric")
  ddHypPart <- DarkDiv(x = dune2, method = "Hypergeometric", removeAbsent = F)
  expect_equal(dimnames(ddHypFull$AllProbs), dimnames(ddHypPart$Pool))
})

