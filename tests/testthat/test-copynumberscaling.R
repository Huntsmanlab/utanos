test_that("Absolute copy number scaling results are correct.", {
  test_data <- readRDS(file = test_path("fixtures", "rascal_test_data.rds"))
  
  rcn <- test_data[["rcn"]]
  
  valid_acn <- test_data[["acn"]]
  valid_best_solns <- test_data[["best_solns"]]
  valid_chosen_solns <- test_data[["chosen_solns"]]
  
  rcn_obj <- DfToQDNAseq(rcn)
  
  solutions <- FindRascalSolutions(rcn_obj)
  
  results <- CalculateACNs(rcn_obj, "mad", solutions, return_sols = TRUE, distance_decimal_places = 7)
  results_acn <- results[["acn_segment_tables"]] %>% dplyr::bind_rows(.id = 'sample_id')
  results_chosen_sols <- results[["rascal_solutions"]]
  
  expect_equal(solutions, valid_best_solns)
  expect_equal(results_chosen_sols, valid_chosen_solns)
  expect_equal(results_acn, valid_acn)
})