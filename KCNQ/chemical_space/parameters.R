

get_system_id <- function(
    user = NULL,
    nodename = NULL,
    verbose = TRUE) {
    if (is.null(user)) {
        user <- Sys.info()[['user']]
    }
    
    if (is.null(nodename)) {
        nodename <- Sys.info()[['nodename']]
    }
    
    system_id <- paste0(user, "@", nodename)
    if (verbose) {
        cat("Using parameters for system '", system_id, "'\n", sep = "")
    }
    system_id
}

get_parameters <- function(
    user = NULL,
    nodename = NULL,
    verbose = FALSE) {

    system_id <- get_system_id(user, nodename, verbose)

    parameters <- list(
        date_code = "20210323",
        project_data_googlesheets_id = "1a7Y93LkvcgQkpqhmQN2CmRNdCjnJiRgIB6IaIaJktKw"
    )

    if (system_id %>% stringr::str_detect("maom@.+[.]arc-ts[.]umich[.]edu")) {
        parameters <- c(
            parameters,
            featurize_substances_program = "~/opt/miniconda3/envs/mplearn/bin/featurize_substances",
            embed_umap_program = "~/opt/miniconda3/envs/mplearn/bin/embed_umap"
        )
    } else if (system_id == "momeara@gimel.cluster.ucsf.bkslab.org") {
        parameters <- c(
            parameters,
            dude_scritps_path = "/mnt/nfs/home/rstein/zzz.scripts/new_DUDE_SCRIPTS/",
            featurize_substances_program = "/nfs/ex9/work/momeara/tools/anaconda3/envs/DeepSEA/bin/featurize_substances",
            embed_umap_program = "/nfs/ex9/work/momeara/tools/anaconda3/envs/DeepSEA/bin/embed_umap")
            
    } else {
        cat("ERROR: Unrecognized system_id='", system_id, "'. Please check parameters.R\n", sep = "")
    }
    parameters
}

parameters <- get_parameters(verbose = TRUE)
