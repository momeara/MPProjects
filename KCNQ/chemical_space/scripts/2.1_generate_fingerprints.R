library(tidyverse)

source("parameters.R")

date_code <- "20220713"


########################
## ECFP4 Fingerprints ##
########################
cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/substances_sanitized_", date_code, ".tsv \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --output_path intermediate_data/project_substances_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --output_path intermediate_data/project_decoys_20210323 \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path /scratch/maom_root/maom99/maom/chembl25_substances.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'substance.smiles' \\
  --output_path intermediate_data/chembl25_substances_20210323 \\
  --verbose
")
cat(cmd, "\n", sep = "")
system(cmd)

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/chembl27_substances_20210323.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'smiles' \\
  --output_path intermediate_data/chembl27_substances_ecfp4_20210323 \\
  --verbose
")
cat(cmd, "\n", sep = "")
system(cmd)



#######################
## APDP Fingerprints ##
#######################
cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/substances_sanitized_", date_code, ".tsv \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/project_substances_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/project_decoys_APDP_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path /scratch/maom_root/maom99/maom/chembl25_substances.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'substance.smiles' \\
  --fingerprint_type 'APDP' \\
  --fingerprint_n_bits 5000 \\
  --output_path intermediate_data/chembl25_substances_APDP_20210323 \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


############################
# Huggingface ChemGTP-4.7M #
############################

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path 'intermediate_data/substances_sanitized_20210323.tsv' \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --fingerprint_type 'huggingface:ncfrey/ChemGPT-4.7M' \\
  --output_path intermediate_data/project_substances_ChemGPT-4.7M_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --fingerprint_type 'huggingface:ncfrey/ChemGPT-4.7M' \\
  --output_path intermediate_data/project_decoys_ChemGPT-4.7M_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/chembl27_substances_20220711.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'smiles' \\
  --fingerprint_type 'huggingface:ncfrey/ChemGPT-4.7M' \\
  --output_path intermediate_data/chembl27_substances_ChemGPT-4.7M_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


############################
# Huggingface ChemGTP-1.2B #
############################

cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path 'intermediate_data/substances_sanitized_20210323.tsv' \\
  --substance_id_field 'substance_name' \\
  --smiles_field 'substance_smiles_rdkit' \\
  --fingerprint_type 'huggingface:/home/maom/opt/huggingface/ChemGPT-1.2B' \\
  --device 'cuda' \\
  --output_path intermediate_data/project_substances_ChemGPT-1.2B_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path intermediate_data/project_decoys_20210323.tsv \\
  --substance_id_field 'substance_zinc_id' \\
  --smiles_field 'substance_smiles' \\
  --fingerprint_type 'huggingface:/home/maom/opt/huggingface/ChemGPT-1.2B' \\
  --device 'cuda' \\
  --output_path intermediate_data/project_decoys_ChemGPT-1.2B_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)


cmd <- paste0(
parameters$featurize_substances_program, " \\
  --library_path raw_data/chembl27_substances_20220711.tsv \\
  --substance_id_field 'zinc_id' \\
  --smiles_field 'smiles' \\
  --fingerprint_type 'huggingface:/home/maom/opt/huggingface/ChemGPT-1.2B' \\
  --device 'cuda' \\
  --output_path intermediate_data/chembl27_substances_ChemGPT-1.2B_", date_code, " \\
  --verbose
")
cat(cmd, sep = "")
system(cmd)

