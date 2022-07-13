


# if the data has already been loaded and then cached in s3, re-download it using these commands
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS12h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS18h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS24h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS3h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS48h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/covid19cq1_SARS_TS6h_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202006_cell_profiler_features/image_scores_CQ1_TS.parquet intermediate_data/")


system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL1_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL1_InfectedCells_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL2_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL2_InfectedCells_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL3_Cell_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/covid19cq1_SARS_TS2PL3_InfectedCells_MasterDataTable.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/image_scores_CQ1_TS_202008.parquet intermediate_data/")
system("aws s3 cp s3://umich-insitro/CQ1/pseudo_time_202008_cell_profiler_features/plate_map_TS_202008.parquet intermediate_data/")
