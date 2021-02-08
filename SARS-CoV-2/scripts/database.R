library(plyr)
library(tidyverse)
library(RMySQL)
library(magrittr)
library(tictoc)


get_test_database_connection <- function(){
    con <- DBI::dbConnect(
        RMySQL::MySQL(),
        host = "covid19cp.cgymeokijgns.us-east-1.rds.amazonaws.com",
        port = 3306,
        user = "covid19cp",
        password = "Genes-brett-Flip-9Bottling")
    con %>% DBI::dbSendQuery("USE covid19cp")
    con
}

set_schema <- function(con, schema) {
    available_schemas <- con %>% DBI::dbGetQuery("SHOW SCHEMAS;")
    if (schema %in% available_schemas$Database) {
        con %>% DBI::dbSendQuery(paste0("USE ", schema))
    } else{
        cat(
            "Unrecognized schema, available schemas are:\n  ",
            available_schemas$Database %>% paste0(collapse = "\n  "), "\n",
            sep = "")
    }
    con
}

# TODO: gather database arguments to be parsed by argparse
# then use these arguments to initialize a database connection
add_mysql_argparse_parameters <- function(parser){
    parser$add_argument(
        "--database_driver",
        action = "store",
        default = "RMySQL",
        dest = "database_driver",
        type = "character")
    parser$add_argument()
        
}

        


get_primary_database_connection <- function(schema="covid19primary") {
    con <- DBI::dbConnect(
        RMySQL::MySQL(),
        host = "covid19primary.cgymeokijgns.us-east-1.rds.amazonaws.com",
        port = 3306,
        user = "covid19primary",
        password = "atop-9most-5Inn-Dandruff9")
#        password = "frightful-bootlace-Cats-flaring8")
    con %>% set_schema(schema)
    con
}


build_indices <- function(plate_id) {
    # if a table is on the LHS of a join
    # an index over the join columns makes the query faster
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX SARS_", plate_id, "_Per_Image_index
       ON ", plate_prefix, "_Per_Image (ImageNumber ASC);"))
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX SARS_", plate_id, "_Per_syncytia_index
       ON ", plate_id, "_Per_syncytia (
            ImageNumber ASC,
            syncytia_Number_Object_Number ASC);"))     
    # 27,299 rows
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX", plate_id, "_Per_Nuclei_index
       ON ", plate_id, "_Per_Nuclei (
            ImageNumber ASC,
            Nuclei_Number_Object_Number ASC);"))
    # <MySQLResult:38546520,0,19>
    
    con %>% DBI::dbSendQuery(paste0("
       CREATE UNIQUE INDEX ", plate_id, "_Per_Cells_index
       ON ", plate_id, "_Per_Cells (
            ImageNumber ASC,
            Cells_Number_Object_Number ASC);"))
    # <MySQLResult:135600456,0,20>
}

#######################


get_acas_database <- function(stage=FALSE){
    con <- DBI::dbConnect(
        RPostgres::Postgres(),
        host=paste0("umsexton", ifelse(stage, "-stage", ""), ".onacaslims.com"),
        port=5432,
        user="acas",
        password="acas")
    con %>% DBI::dbSendQuery("SET search_path TO acas;")
    con
}


