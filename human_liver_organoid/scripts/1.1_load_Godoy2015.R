library(tidyverse)
library(fs)
library(monocle3)
library(MPStats)
library(Seurat)

source("parameters.R")

################3
# try to read the bigwig files directly from GEO
library("rtracklayer")
data <- data.frame(
    bw_path = Sys.glob(
        paste0(parameters$data_base_dir, "/scRNAseq/gao2017/*.bw"))) %>%
    dplyr::mutate(
        gsm_id = bw_path %>% stringr::str_extract("GSM[0-9]+"),
        run_id = bw_path %>% stringr::str_extract("GSM.*$") %>%
            stringr::str_replace("^GSM[0-9]+_", "") %>%
            stringr::str_replace("_mRNA[.]bw$", ""))

# this gives a GRanges file but I have no idea what gene model was used?
z <- rtracklayer::import.bw(data$bw_path[1])

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txbygene <- transcriptsBy(txdb, "gene")
overlaps = findOverlaps(z, txbygene)

txdb <- GenomicFeatures::makeTxDbFromBiomart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl")


####################################################
# try to re-process the sra files with RSEM

sra_dir <- paste0(parameters$scratch_dir, "/sra")
if (!dir.exists(sra_dir)) {
    dir.create(sra_dir, recursive = TRUE)
}

system(paste0("cd ", sra_dir, " && prefetch SRR5974287 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974288 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974289 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974290 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974291 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974292 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974293 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974294 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974295 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974296 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974297 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974298 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR5974299 --progress 1"))
system(paste0("cd ", sra_dir, " && prefetch SRR6179879 --progress 1"))

system(paste0(
parameters$rsem_prepare_reference_program, " \\
  --num-threads 4 \\
  --bowtie2 --bowtie2-path ", parameters$bowtie2_path, " \\
  ", parameters$scratch_dir, "/ref_genome_GRCh38/gencode.v38.transcripts.fa \\
  ", parameters$scratch_dir, "/ref_genome_GRCh38
"))

tag <- "gao2017"
n_jobs <- 6
job_dir <- paste0(parameters$scratch_dir, "/estimated_expression_", tag)
if (!dir.exists(job_dir)) {
    cat("Creating job directory: ", job_dir, "\n", sep = "")
    dir.create(job_dir, recursive = TRUE)
}
base_dir <- getwd()

cmd_str <- paste0(
    "sbatch ",
    "--account=", parameters$slurm_account, " ",
    "--mail-user=", parameters$slurm_mail_user, " ",
    "--mail-type=BEGIN,END,FAIL ",
    "--array=1-", n_jobs, " ",
    "--output='", job_dir, "/%j.log' ",
    "--time=12:00:00 ",
    "--nodes=", n_jobs, " ",
    "--export=",
    "TAG='", tag, "',",
    "BASE_DIR='", base_dir, "',",
    "JOB_DIR='", job_dir, "' ",
    "scripts/run_estimate_expression_SLURM_wrapper.sh")
info_message <- "Monitor progress with 'squeue'"

cat(cmd_str, "\n")
system(cmd_str)
cat(info_message, "\n", sep = "")
cat("Check results when done: intermediate_data/estimated_expression_", tag, "/logs\n", sep = "")

####################################

samples_gao2017 <- readr::read_tsv(
    file = paste0(
        parameters$data_base_dir,
        "/scRNAseq/gao2017/estimated_expression_gao2017/samples.tsv"))

library(tximport)
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#RSEM

data_gao2017 <- tximport::tximport(
    files = paste0(
        parameters$data_base_dir,
        "/scRNAseq/gao2017/estimated_expression_gao2017/",
        c(
        "SRR5974287.genes.results",
        "SRR5974288.genes.results",
        "SRR5974289.genes.results",
        "SRR5974290.genes.results",
        "SRR5974291.genes.results",
        "SRR5974292.genes.results",
        "SRR5974293.genes.results",
        "SRR5974294.genes.results",
        "SRR5974295.genes.results",
        "SRR5974296.genes.results",
        "SRR5974297.genes.results",
        "SRR5974298.genes.results",
        "SRR5974299.genes.results",
        "SRR6179879.genes.results")),
    txIn = FALSE,
    type = "rsem")

transcripts_gao2017 <- data.frame(
    hg38_transcript_id = rownames(data_gao2017$abundance)) %>%
    tidyr::separate(
        col = hg38_transcript_id,
        into = c(
            "ensembl_transcript",
            "ensembl_gene_version",
            "havana_gene",
            "havana_transcript",
            "transcript_start",
            "gene_name",
            "transcript_length",
            "transcript_type",
            "extra"),
        sep = "[|]",
        remove = FALSE) %>%
    dplyr::mutate(
        ensembl_gene = ensembl_gene_version %>%
            stringr::str_extract("^ENSG[0-9]+"))

genes_gao2017 <- transcripts_gao2017 %>%
    dplyr::distinct(ensembl_gene)


# off chip
load(file = "intermediate_data/cds_CZ.Rdata")
# on chip
load(file = "intermediate_data/cds_CZ_2.Rdata")

cds <- monocle3::combine_cds(
    cds_list = list(
        off_chip = cds_CZ,
        on_chip = cds_CZ_2))

cds <- cds %>%
    MPStats::compute_qc_covariates()

cds <- cds[
    # genes that are observed in at least 5 cells
    SingleCellExperiment::counts(cds) %>% Matrix::rowSums() >= 5,
    # cells that have at least 1000 reads which are at most 30% mitochondrial
    (SummarizedExperiment::colData(cds)[["count_depth"]] >= 10000) &
    (SummarizedExperiment::colData(cds)[["mt_fraction"]] < .3)]

# # get warning and errors about the names of the reducedDim objects (e.g. PCA, UMAP) names
# # so remote them for now
# reducedDims(cds) <- NULL
# # get error about missing logcounts, resolve by setting data = NULL
# # https://github.com/satijalab/seurat/issues/3746#issuecomment-731419868
# cds <- Seurat::as.Seurat(cds, data = NULL)
#
# Seurat::Idents(cds) <- cds@meta.data$sample
#
# cds <- cds %>% Seurat::SCTransform(
#     assay = "originalexp",
#     vst.flavor = "v2")
#

genes_CZ <- data.frame(
    ensembl_gene = cds %>%
        SingleCellExperiment::counts() %>%
        rownames())

genes_all <- genes_gao2017 %>%
    dplyr::mutate(in_gao2017 = TRUE) %>%
    dplyr::full_join(
        genes_CZ %>%
        dplyr::mutate(in_CZ = TRUE),
        by = "ensembl_gene")
#   in_gao2017 in_CZ     n
# 1       TRUE  TRUE 23924
# 2       TRUE    NA 36681
# 3         NA  TRUE    62

hsapiens_genes <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
    mart = biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"))

genes_common <- genes_gao2017 %>%
    dplyr::inner_join(genes_CZ, by = "ensembl_gene") %>%
    dplyr::left_join(
        hsapiens_genes,
        by = c("ensembl_gene" = "ensembl_gene_id"))


agg_gao2017 <- genes_common %>%
    dplyr::inner_join(
        dplyr::bind_cols(
            transcripts_gao2017 %>%
            dplyr::select(ensembl_gene),
            data_gao2017$counts %>% data.frame()) %>%
        dplyr::group_by(ensembl_gene) %>%
        dplyr::summarize(
            dplyr::across(dplyr::everything(), ~sum(.x))),
        by = "ensembl_gene")
names(agg_gao2017) <- c(
    "ensembl_gene", "hgnc_symbol", "description",
    paste0(samples_gao2017$sample, "_", samples_gao2017$replicate))
agg_gao2017 %>%
    dplyr::left_join(genes_common, by = "ensembl_gene") %>%
    readr::write_tsv(
        "product/figures/CZ_CZ_2/gao2017_20210720.tsv")



agg_CZ <- genes_common %>%
    dplyr::left_join(
        dplyr::bind_cols(
            genes_CZ,
            cds %>% SingleCellExperiment::counts() %>% as.data.frame()),
        by = "ensembl_gene")

n_pseudo_bulk_samples <- 10

CZ_sample_ids <- data.frame(
    sample_id = cds %>% SingleCellExperiment::counts() %>% colnames()) %>%
    dplyr::mutate(
        barcode = sample_id %>% stringr::str_extract("^[ACGT]+[-][0-9]+"),
        sample = sample_id %>% stringr::str_replace("^[ACGT]+[-][0-9]+_", "")) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
        pseudo_bulk_id = 1:n_pseudo_bulk_samples %>%
            rep(length.out = dplyr::n()) %>%
            sample(size = dplyr::n())) %>%
    dplyr::ungroup()

pb_CZ <- dplyr::bind_cols(
    CZ_sample_ids,
    agg_CZ %>%
    dplyr::select(-ensembl_gene, -hgnc_symbol, -description) %>%
    t() %>% data.frame) %>%
    dplyr::select(-sample_id, -barcode) %>%
    dplyr::group_by(sample, pseudo_bulk_id) %>%
    dplyr::summarize(dplyr::across(everything(), ~sum(.x))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_id = paste0(sample, "_", pseudo_bulk_id)) %>%
    tibble::column_to_rownames("sample_id") %>%
    dplyr::select(-sample, -pseudo_bulk_id) %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(
        ensembl_gene = agg_CZ$ensembl_gene,
        hgnc_symbol = agg_CZ$hgnc_symbol,
        description = agg_CZ$description,
        .before = 1)

pb_CZ %>%
    dplyr::left_join(genes_common, by = "ensembl_gene") %>%
    readr::write_tsv(
        "product/figures/CZ_CZ_2/on_off_chip_pseudobulk_20210720.tsv")

countData <- dplyr::bind_cols(
    agg_gao2017 %>% dplyr::select(-ensembl_gene, -hgnc_symbol, -description),
    pb_CZ %>% dplyr::select(-ensembl_gene, -hgnc_symbol, -description)) %>%
    as.matrix(rownames.value = genes_common$ensembl_gene) %>%
    round()
rownames(countData) <- genes_common$ensembl_gene

colData <- dplyr::bind_rows(
    samples_gao2017 %>% dplyr::select(
        sample,
        replicate),
    data.frame(
        sample = names(pb_CZ)[2:21] %>% stringr::str_replace("_[0-9]+$", ""),
        replicate = names(pb_CZ)[2:21] %>% stringr::str_extract("[0-9]+$") %>% as.numeric()))

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData = colData,
    design = ~ sample)

dds <- DESeq2::DESeq(dds)

contrasts <- dplyr::bind_rows(
    data.frame(numerator = "on_chip", denominator = "off_chip"),
    data.frame(numerator = "off_chip", denominator = "hiHep"),
    data.frame(numerator = "off_chip", denominator = "hiPSC"),
    data.frame(numerator = "off_chip", denominator = "PHH"),
    data.frame(numerator = "off_chip", denominator = "UCF"),
    data.frame(numerator = "on_chip", denominator = "hiHep"),
    data.frame(numerator = "on_chip", denominator = "hiPSC"),
    data.frame(numerator = "on_chip", denominator = "PHH"),
    data.frame(numerator = "on_chip", denominator = "UCF"),
    data.frame(numerator = "PHH", denominator = "hiHep"),
    data.frame(numerator = "PHH", denominator = "hiPSC"),
    data.frame(numerator = "PHH", denominator = "UCF"),
    data.frame(numerator = "hiHep", denominator = "hiPSC"),
    data.frame(numerator = "hiPSC", denominator = "UCF"),
    data.frame(numerator = "hiHep", denominator = "UCF")) %>%
    dplyr::mutate(condition = "sample")

results <- contrasts %>%
    dplyr::rowwise() %>%
    dplyr::do({
        contrast <- .
        results <- DESeq2::results(
            object = dds,
            contrast = c(contrast$condition, contrast$numerator, contrast$denominator))
        ensembl_gene <- rownames(results)
        results %>%
            data.frame() %>%
            dplyr::mutate(
                numerator = contrast$numerator,
                denominator = contrast$denominator,
                ensembl_gene = ensembl_gene) %>%
            dplyr::left_join(
                genes_common,
                by = "ensembl_gene")
    }) %>%
    dplyr::ungroup()


results %>%
    readr::write_tsv("product/figures/CZ_CZ_2/compare_gao2017/differential_expression_on_off_chip_pseudobulk_gao2017_20210720.tsv")

results_count <- results %>%
    dplyr::filter(abs(log2FoldChange) > 2, padj < 1e-20) %>%
    dplyr::count(numerator, denominator)

# ####
# TXNIP
# TRIB3
# DDIT3
# EGFL6
# NNMT
# RSPO2
# TAGLN
# DSC1

marker_genes <- readr::read_tsv("raw_data/marker_genes_20210306.tsv")

liver_markers <- results %>%
    dplyr::filter(baseMean > 100) %>%
    dplyr::inner_join(
        marker_genes,
        by = c("hgnc_symbol" = "gene"))

liver_markers %>%
    readr::write_tsv("product/figures/CZ_CZ_2/compare_gao2017/de_liver_markers_20210726.tsv")

liver_markers %>%
    dplyr::group_by(Source) %>%
    dplyr::do({
        markers <- .
        marker_source <- markers$Source[1]
        cat("Making plot for ", marker_source, "\n", sep = "")

        de_not_significant <- results %>%
            dplyr::filter(baseMean > 100) %>%
            dplyr::group_by(numerator, denominator) %>%
            dplyr::filter(
                (abs(log2FoldChange) <= 1) |
                (log10(padj) > -5)) %>%
            dplyr::ungroup()

        de_significant <- results %>%
            dplyr::filter(baseMean > 100) %>%
            dplyr::group_by(numerator, denominator) %>%
            dplyr::filter(
                (abs(log2FoldChange) > 1) &
                (log10(padj) < -5)) %>%
            dplyr::ungroup()

        ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::geom_point(
                data = de_not_significant,
                mapping = ggplot2::aes(
                    x = log2FoldChange,
                    y = -log10(padj)),
                size = .4,
                color = "grey50",
                alpha = .3) +
            ggplot2::geom_point(
                data = de_significant,
                mapping = ggplot2::aes(
                    x = log2FoldChange,
                    y = -log10(padj)),
                color = "darkgreen",
                size = .4,
                alpha = .5) +
            ggplot2::geom_point(
                data = markers,
                mapping = ggplot2::aes(
                    x = log2FoldChange,
                    y = -log10(padj),
                    color = gene_set),
                size = 1,
                alpha = 1) +
            ggrepel::geom_text_repel(
                data = markers,
                mapping = ggplot2::aes(
                    x = log2FoldChange,
                    y = -log10(padj),
                    color = gene_set,
                    label = hgnc_symbol),
                size = 1.8,
                force = 10,
                min.segment.length = 0) +
            ggplot2::ggtitle(
                label = "Compare Gao2017 and on- and off-chip pseudo bulk expression with DESeq2",
                subtitle = paste0("Marker set: ", marker_source)) +
            ggplot2::scale_x_continuous(expression("Average log[2] FC")) +
            ggplot2::scale_y_continuous(
                expression("-log[10] Corrected P-value")) +
            ggplot2::facet_wrap(facets = dplyr::vars(numerator, denominator)) +
            ggplot2::theme(legend.position = "bottom")

        ggplot2::ggsave(
            filename = paste0(
                paste0("product/figures/CZ_CZ_2/compare_gao2017/de_volcano_", fs::path_sanitize(marker_source), "_20210726.pdf")),
            width = 15,
            heigh = 15,
            useDingbats = FALSE)
    })
