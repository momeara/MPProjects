

parameters <- list(
    base_dir = ".",

    # this is the Sexton Lab Turbo drive as it is mounted on the UMich great-lakes cluster
    # link valid as of 2/8/21
    data_base_dir = "/nfs/turbo/umms-jzsexton",

    # HLOs/datasets
    # a table of the datasets that have been collected for the project
    datasets_googlesheet_id = "1R_dd-ccWajPZxOWkcjfX8EWNkAalNHmq0HfWn75-8Hk",
    datasets_googlesheet_name = "Summary",

    # path where the datasets table is stored
    datasets_fname = "raw_data/datasets_20210208.tsv",

    marker_genes_googlesheet_id = "1R_dd-ccWajPZxOWkcjfX8EWNkAalNHmq0HfWn75-8Hk",
    marker_genes_googlesheet_sheet = "Marker Genes",

    marker_genes_fname = "raw_data/marker_genes_20210306.tsv",

    scratch_dir = "/scratch/maom_root/maom99/maom/HLOs",

    # cluster scratch used for computing expression profiles
    expression_scratch_dir = "/scratch/maom_root/maom99/maom",
    reference_genome_path = "/scratch/maom_root/maom99/maom/HLOs/ref_genome_GRCh38",
    bowtie2_path = "/home/maom/opt/bowtie2-2.3.4.1-linux-x86_64",
    fastq_dump_program = "/home/maom/opt/bin/fastq-dump",
    rsem_prepare_reference_program = "/home/maom/opt/RSEM/rsem-prepare-reference",
    rsem_calculate_expression_program = "/home/maom/opt/RSEM/rsem-calculate-expression",
    slurm_account = "maom99",
    slurm_mail_user = "maom@umich.edu",
    slurm_partition = "standard"
)
