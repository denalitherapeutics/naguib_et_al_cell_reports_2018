#' Retrieve raw abundance data from NCBI GEO
#'
#' @note An internet connection is required to retrieve the processed counts
#' from NCBI's GEO repository
#' @param dataset Scalar character, the experiment of interest, either `HEK293`
#' or `Fibroblasts`.
#' @importFrom edgeR DGEList
#' @importFrom dplyr bind_rows
#' @return [edgeR::DGEList] object
#' @export
load_dataset <- function(
  experiment = c("HEK293", "Fibroblasts", "HEK293_71D")) {
  gene_anno_columns <- c("entrez_id", "hgnc_symbol", "gene_biotype",
                         "ensembl_gene_id")

  experiment <- match.arg(experiment)
  accession <- switch(
    experiment,
    "HEK293" = "GSE116267",
    "Fibroblasts" = "GSE116266",
    "HEK293_71D" = "GSE121757"
  )
  count_file <- switch(
    experiment,
    "HEK293" = system.file(package = "SUPT4H1", "extdata",
                           "GSE116267_SUPT4H1_HEK293.tab.gz", mustWork = TRUE),
    "Fibroblasts" = system.file(package = "SUPT4H1", "extdata",
                                "GSE116266_SUPT4H1_Fibroblasts.tab.gz",
                                mustWork = TRUE),
    "HEK293_71D" = system.file(package = "SUPT4H1", "extdata",
                               "GSE121757_SUPT4H1_HEK293_71D.tab.gz",
                               mustWork = TRUE)
  )
  counts <- read.delim(count_file, row.names = 1, stringsAsFactors = FALSE,
                       check.names = FALSE)

  sample_anno_file <- switch(
    experiment,
    "HEK293" = system.file(package = "SUPT4H1", "extdata",
                           "GSE116267_SUPT4H1_HEK293.csv", mustWork = TRUE),
    "Fibroblasts" = system.file(package = "SUPT4H1", "extdata",
                                "GSE116266_SUPT4H1_Fibroblasts.csv",
                                mustWork = TRUE),
    "HEK293_71D" = system.file(package = "SUPT4H1", "extdata",
                                "GSE121757_SUPT4H1_HEK293_71D.csv",
                                mustWork = TRUE)
  )
  sample_anno <- read.csv(sample_anno_file, stringsAsFactors = FALSE,
                          row.names = "sample_id")

  # split off the gene annotations
  gene_anno <- counts[, gene_anno_columns]
  gene_anno$spikein <- grepl("^ERCC-", row.names(gene_anno))

  # annotate ERCC spike-in controls
  spikes <- ercc()
  spikes$ercc_logFC <- switch(
    experiment,
    "HEK293" = round(log2(spikes$conc_mix2 / spikes$conc_mix1), digits = 2),
    "Fibroblasts" = round(
      log2(spikes$conc_mix2 / spikes$conc_mix1), digits = 2),
    "HEK293_71D" = round(log2(spikes$conc_mix2 / spikes$conc_mix1), digits = 2)
  )
  spikes <- merge(
    x = gene_anno[gene_anno$spikein, ],
    y = spikes[, c("ercc_id", "conc_mix1", "conc_mix2", "ercc_logFC")],
    by.x = "ensembl_gene_id", by.y = "ercc_id", all.x = TRUE)
  gene_anno <- dplyr::bind_rows(gene_anno[!gene_anno$spikein, ], spikes)
  row.names(gene_anno) <- gene_anno$ensembl_gene_id

  # create DGEList object
  counts <- data.matrix(counts[, !colnames(counts) %in% gene_anno_columns])
  dge <- edgeR::DGEList(
    counts = counts,
    samples = sample_anno[colnames(counts), ],
    genes = gene_anno[row.names(counts), ]
  )
  if (experiment == "Fibroblasts") {
    dge$samples$replicate <- paste0("R", dge$samples$replicate)
  }
  if (experiment == "HEK293") {
    dge$samples$sample_id <- colnames(dge)
  }
  return(dge)
}

#' Load precomputed table of intronic counts for HEK293 cells
#'
#' @importFrom edgeR DGEList calcNormFactors
#' @return [edgeR::DGEList] object
#' @export
load_intron_counts <- function() {
  gene_anno_columns <- c("entrez_id", "hgnc_symbol", "gene_biotype",
                         "ensembl_gene_id")

    count_file <- system.file(
    "extdata", "GSE116267_SUPT4H1_HEK293_introns.csv.gz", package = "SUPT4H1",
    mustWork = TRUE)
  counts <- read.csv(count_file)
  row.names(counts) <- counts$ensembl_gene_id

  sample_anno_file <- system.file(package = "SUPT4H1", "extdata",
                           "GSE116267_SUPT4H1_HEK293.csv", mustWork = TRUE)
  sample_anno <- read.csv(sample_anno_file, stringsAsFactors = FALSE,
                          row.names = "sample_id")

  # split off the gene annotations
  gene_anno <- counts[, gene_anno_columns]

  # create DGEList object
  counts <- data.matrix(counts[, !colnames(counts) %in% gene_anno_columns])
  dge <- edgeR::DGEList(
    counts = counts,
    samples = sample_anno[colnames(counts), ],
    genes = gene_anno[row.names(counts), ]
  )
  dge <- edgeR::calcNormFactors(dge)
  return(dge)
}

#' Load flow cytometry results
#'
#' @param experiment Scalar character, one of `RNAi` or `RNAse`, specifying the
#' results from cells treated with SUPT4H1 / control siRNAs or RNAse,
#' respectively.
#' @return data.frame with columns `signal` and `condition`.
#' @export
#' @examples
#' head(load_flow_data("siRNA"))
load_flow_data <- function(experiment = c("siRNA", "RNAse")) {
  experiment <- match.arg(experiment)
  csv_file <- switch(
    experiment,
    "siRNA" = "HEK293_SUPT4H1_RNAi.csv.gz",
    "RNAse" = "HEK293_FITC_RNAse.csv.gz"
  )
  csv_path <- system.file("extdata", csv_file, package = "SUPT4H1",
                          mustWork = TRUE)
  read.csv(csv_path, stringsAsFactors = FALSE)
}
