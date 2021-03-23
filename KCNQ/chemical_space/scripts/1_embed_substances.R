
system("
date_code=$(date '+%Y%m%d')
embed_substances \
        --library_path ../docking/databases/project_20201214/ \
        --output_path intermediate_data/project_substances_${date_code}
")
