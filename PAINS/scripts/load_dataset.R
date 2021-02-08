
load_dataset <- function (
    dataset_id,
    verbose = FALSE) {
    if (verbose) {
        cat("Loading dataset '", dataset_id, "'\n", sep = "")
    }

    # setup local paths
    dataset_path <- paste0("raw_data/", dataset_id)
    if (dir.exists(dataset_id)) {
        cat("WARNING: '", dataset_path, "' already exists\n", sep = "")
        command <- paste0("rm -rf ", dataset_path, sep = "")
        if (verbose) {
            cat(command, "\n", sep = "")
        }
        system(command)
    } else {
        command <- paste0("mkdir ", dataset_path)
        if (verbose) {
            cat(command, "\n", sep = "")
        }
        system(command)
    }

    command <- paste0("sudo cp \\
	 ~/bucket_cellprofilerdata/PAINS/", dataset_id, "/cpdata.h5 \\
	 ", dataset_path, "/cpdata.h5
    sudo chmod a+r ", dataset_path, "/cpdata.h5")
    if (verbose) {
        cat(command, "\n", sep = "")
    }
    system(command)
}
