
system("
date_code=$(date '+%Y%m%d')
embed_substances \
        --library_path ../data_curation/substances.tsv/ \
        --output_path intermediate_data/project_substances_${date_code}
")
