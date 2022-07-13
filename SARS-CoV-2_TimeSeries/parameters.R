

parameters <- list(
    databases = list(
        covid19cq1 = list(
            database_type = "mysql",
            host = "covid19primary.cgymeokijgns.us-east-1.rds.amazonaws.com",
            port = 3306,
            user = "covid19primary",
            password = "atop-9most-5Inn-Dandruff9",
            schema = "covid19cq1"
        )
    ),
    python_env = "/home/ubuntu/anaconda3/envs/sextonlab",

    base_dir = "/home/ubuntu/projects/SARS-CoV-2_TimeSeries",

    plate_ids_fname = "raw_data/plate_ids_20210125.tsv",

    # google drive path:
    # momeara_lab/Morphological Profiling/SARS-CoV-2/Screening Runs/SARS Project Plate Log
    plate_map_googlesheet_id = "1sJ7UekU9aEh-56Sy4mpeOmVCGlEh8wcnKjj6ubtBI-o"
    )

