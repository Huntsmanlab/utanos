test_that("Segments are properly converted from relative to absolute positions.", {
  rel_seg_pos_all_chr <- readRDS(file = "./tests/testthat/fixtures/rel_seg_pos_all_chr.rds")

  chromosomes <- rel_seg_pos_all_chr$chromosomes
  rel_start_pos <- rel_seg_pos_all_chr$rel_start_pos
  rel_end_pos <- rel_seg_pos_all_chr$rel_end_pos
  build = "GRCh37"

  expect_equal(RelToAbsSegPos(chromosomes = chromosomes, rel_start_pos = rel_start_pos, rel_end_pos = rel_end_pos, build = build), rel_seg_pos_all_chr$result)
})
