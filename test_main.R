#!/usr/bin/Rscript
source("main.R")
library(testthat)

test_that("loading csv works using load_expression()", {
  result_tib <- load_expression(unzip("data/example_intensity_data.zip"))
  expect_equal(dim(result_tib), c(54675, 36))
  expect_true(is_tibble(result_tib))
})

test_that("the correct rows are filtered out in filter_15()", {
  test_tib <- tibble(probeids=c('1_s_at', '2_s_at', '3_s_at', '4_s_at'),
                     GSM1=c(1.0, 3.95, 4.05, 0.5),
                     GSM2=rep(1.6, 4),
                     GSM3=rep(2.5, 4),
                     GSM4=rep(3.99, 4),
                     GSM5=rep(3.0, 4),
                     GSM6=rep(1.0, 4),
                     GSM7=rep(0.5, 4))
  expect_equal(pull(filter_15(test_tib)), c("2_s_at", "3_s_at"))
})

test_that("affy ids can be converted to HGNC names properly using affy_to_hgnc()", {
  # biomaRt super buggy so we can try to account for not connecting well
  response <- try(affy_to_hgnc(tibble('1553551_s_at')), TRUE)
  if (grepl("Error", response[1])) {
    expect_warning(warning("Could not connect to ENSEMBL."))
  } 
  else {
    expect_equal(response$hgnc_symbol, c("MT-ND1", "MT-TI", "MT-TM", "MT-ND2", "MT-TW"))
  }
})

test_that("reduce_data() is correctly changing the size and shape of the tibble", {
  t_tibble <- tibble(probeids = c("1_s_at", "2_s_at", "3_s_at"),
                     GSM1 = c(9.5, 7.6, 5.5),
                     GSM2 = c(9.7, 7.2, 2.9),
                     GSM3 = c(6.9, 4.3, 6.8))
  names <- tibble(affy_hg_u133_plus_2 = c("1_s_at", "2_s_at", "3_s_at"),
                  hgnc_symbol = c("A-REAL-GENE", "SONIC", "UTWPU"))
  good <- c("A-REAL-GENE")
  bad <- c("SONIC")
  reduce_test <- reduce_data(t_tibble, names, good, bad)
  result <- tibble(probeids = c("1_s_at", "2_s_at"),
                   hgnc_symbol = c("A-REAL-GENE", "SONIC"),
                   gene_set = c("good", "bad"),
                   GSM1 = c(9.5, 7.6),
                   GSM2 = c(9.7, 7.2),
                   GSM3 = c(6.9, 4.3))
  expect_equal(reduce_test, result)
})

test_that("plot_ggplot() correctly creates a boxplot from sample data", {
  plot_tib <- tibble(probeids = c("202274_at", "202541_at", "202542_s_at", "203919_at"),
                     hgnc_symbol = c("ACTG2", "AIMP1", "AIMP1", "TCEA2"),
                     gene_set = rep("good", 4),
                     GSM1 = c(8.05, 8.40, 9.55, 4.44),
                     GSM2 = c(7.74, 7.11, 8.48, 5.39))
  p <- plot_ggplot(plot_tib)
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomBoxplot")
})
