#!/bin/bash

rm -rf decoys
rm -rf decoys_final


cp ../project_20210323/substances.smi .
SUBSTANCES_FNAME=substances.smi

time python ../../scripts/dude_scripts/0000_protonate_setup_dirs.py ${SUBSTANCES_FNAME} decoys
python ../../scripts/dude_scripts/0001_qsub_generate_decoys.py decoys
python ../../scripts/dude_scripts/0002_qsub_filter_decoys.py decoys
time python ../../scripts/dude_scripts/0003_copy_decoys_to_new_dir.py decoys decoys_final


bash find $(pwd)/decoys_final/decoys -name '*db2.gz' > database.sdi
