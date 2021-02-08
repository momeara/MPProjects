

parameters <- list(
    databases = list(
        covid18cq1 = list(
            host = "covid19primary.cgymeokijgns.us-east-1.rds.amazonaws.com",
            port = 3306,
            user = "covid19primary",
            password = "atop-9most-5Inn-Dandruff9",
            schema = "covid19cq1"
        ),

    python_env = "/home/ubuntu/anaconda3/envs/sextonlab",
        
    base_dir = "/home/ubuntu/projects/SARS-CoV-2_TimeSeries",

    plate_ids_fname = "raw_data/plate_ids_20210125.tsv",
    
    # this is a downloaded copy of the Google drive "SARS Project Plate Log"
    plate_map_fname = "raw_data/plate_map_20210125.xlsx"


    )

