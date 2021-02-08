

run_analysis:
	Rscript scripts/1.1_read_well_scores.R
	Rscript scripts/1.2_read_compound_moa.R
	Rscript scripts/2_model_cell_count_by_batch_vars.R
	Rscript scripts/2_model_score_by_batch_vars.R
	Rscript scripts/2_model_score_by_cell_count.R
	Rscript scripts/3_drc_dose_response.R

clean:
	rm -rf intermediate_data
	mkdir intermediate_data
	rm -rf product
	mkdir product
