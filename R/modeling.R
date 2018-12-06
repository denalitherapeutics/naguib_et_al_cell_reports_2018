#' Fit a mixed linear model to RNA-seq data from HEK293 cells
#'
#' Differential expression analysis using functions and methods from the
#' `edgeR` and `limma` packages.
#'
#' 1. The [edgeR::filterByExpr()] method is called to determine which genes have
#'    sufficiently large counts to be retained in a statistal analysis, based on
#'    the experimental design.
#' 2. The [edgeR::calcNormFactors()] method is used to calculate a TMM scaling
#'    factor, either based on all endogenous genes or only on the external ERCC
#'    spike-ins retained after filtering, as determined by the `normalize`
#'    argument.
#' 3. The [limma::voom()] function is used to transform raw count data to
#'    log2-counts per million (logCPM), estimate the mean-variance
#'    relationship and compute appropriate observation-level weights.
#' 4. The [limma::duplicateCorrelation()] function is used to estimate the
#'    correlation between between technical replicates (e.g. within each
#'    cell line).
#' 5. The [limma::voom()] function is called a second time, this time also
#'    providing the average estimated inter-duplicate correlation determined
#'    above.
#' 6. A linear model is fit for each gene with the [limma::lmFit()] function.
#' 7. Moderated t-statistics and log-odds of differential expression are
#'    calculated with the [limma::eBayes()] function using an empirical Bayes
#'    approach.
#' A linear mixed model is fit by REML individually for each gene, taking
#' account the correlation between technical replicates performed within each
#' cell line.
#'
#' @param dge DGEList with raw counts
#' @param normalize Scalar string, calculate a TMM scaling factor using either
#' `spike-in` or `global` normalization.
#' @param output Scalar character specifying which output to return, either `v`
#' to return output of [limma::voom()] or `fit` to return the output of
#' [limma::eBayes()].
#' @importFrom limma voom duplicateCorrelation lmFit eBayes
#' @importFrom edgeR filterByExpr calcNormFactors [.DGEList
#' @return [limma::MArrayLM] fit object
#' @export
#' @examples
#' library(limma)
#' # global shift is detected with spike-in normalization
#' fit_spikein <- fit_hek293(normalize = "spike-in")
#'
#' limma::plotMA(fit_spikein, coef = "SUPT4H1",
#'     main = "HEK293: SUPT4H1 RNAi vs Scramble (spike-in normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
#' # but not with global normalization
#' fit_global <- fit_hek293(normalize = "global")
#' limma::plotMA(fit_global, coef = "SUPT4H1",
#'     main = "HEK293: SUPT4H1 RNAi vs Scramble (global normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
fit_hek293 <- function(dge, normalize = c("spike-in", "global", "none"),
                       output = c("fit", "v")) {
  normalize <- match.arg(normalize)
  output <- match.arg(output)

  if (missing(dge)) {
    dge <- load_dataset(experiment = "HEK293")
  }

  # design matrix
  design <- model.matrix(~ treatment + cell_line, data = dge$samples)
  colnames(design) <- gsub("^treatment", "", colnames(design))

  # filtering
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, ]

  # normalization
  dge$genes$spikein <- grepl("^ERCC-", row.names(dge))
  if (normalize == "global") {
    dge_global <- edgeR::calcNormFactors(dge[!dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_global$samples$norm.factors
  } else if (normalize == "spike-in") {
    dge_spikein <- edgeR::calcNormFactors(dge[dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_spikein$samples$norm.factors
  } else {
    message("Using size factors as they were provided in the DGEList input!")
  }

  # voom transformation (without correlation information)
  v <- limma::voom(
    counts = dge,
    design = design,
    plot = FALSE,
    save.plot = FALSE
  )
  # Model duplicate correlation on previously voom-transformed data
  blocks <- with(v$targets, paste(cell_line, treatment))
  corfit <- suppressWarnings(
    limma::duplicateCorrelation(v[v$genes$spikein == FALSE, ],
                                design, block = blocks))

  # secound round of voom transformation, now with correlation information
  v <- limma::voom(
    counts = dge,
    design = design,
    correlation = corfit$consensus,
    plot = FALSE,
    save.plot = FALSE
  )
  if (output == "v") {
    return(v)
  }

  # fit linear model
  fit <- limma::lmFit(object = v[, row.names(design)], design = design,
                      block = blocks, correlation = corfit$consensus)
  fit <- limma::eBayes(fit)
  return(fit)
}

#' Fit a linear model to RNA-seq data from a single HEK293 reporter cell line
#'
#' Differential expression analysis using functions and methods from the
#' `edgeR` and `limma` packages.
#'
#' 1. The [edgeR::filterByExpr()] method is called to determine which genes have
#'    sufficiently large counts to be retained in a statistal analysis, based on
#'    the experimental design.
#' 2. The [edgeR::calcNormFactors()] method is used to calculate a TMM scaling
#'    factor, either based on all endogenous genes or only on the external ERCC
#'    spike-ins retained after filtering, as determined by the `normalize`
#'    argument.
#' 3. The [limma::voom()] function is used to transform raw count data to
#'    log2-counts per million (logCPM), estimate the mean-variance
#'    relationship.
#' 4. A linear model is fit for each gene with the [limma::lmFit()] function.
#' 5. Moderated t-statistics and log-odds of differential expression are
#'    calculated with the [limma::eBayes()] function using an empirical Bayes
#'    approach.
#' A linear  model is fitindividually for each gene.
#'
#' @param dge DGEList with raw counts
#' @param line Scalar string, the name of the cell line to model
#' @param normalize Scalar string, calculate a TMM scaling factor using either
#' `spike-in` or `global` normalization. If `none`, the `dge` object is expected
#' to be normalized already.
#' @param output Scalar character specifying which output to return, either `v`
#' to return output of [limma::voom()] or `fit` to return the output of
#' [limma::eBayes()].
#' @importFrom limma voom duplicateCorrelation lmFit eBayes
#' @importFrom edgeR filterByExpr calcNormFactors [.DGEList
#' @return [limma::MArrayLM] fit object
#' @export
#' @examples
#' library(limma)
#' # global shift is detected with spike-in normalization
#' fit_spikein <- fit_hek293_lines(normalize = "spike-in", line = "Control")
#'
#' limma::plotMA(fit_spikein, coef = "SUPT4H1",
#'     main = "HEK293: SUPT4H1 RNAi vs Scramble (spike-in normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
#' # but not with global normalization
#' fit_global <- fit_hek293(normalize = "global")
#' limma::plotMA(fit_global, coef = "SUPT4H1",
#'     main = "HEK293: SUPT4H1 RNAi vs Scramble (global normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
fit_hek293_lines <- function(dge, line = c("Control", "DARP1", "DARP2"),
                             normalize = c("spike-in", "global"),
                             output = c("fit", "v")) {
  line <- match.arg(line)
  normalize <- match.arg(normalize)
  output <- match.arg(output)

  if (missing(dge)) {
    dge <- load_dataset(experiment = "HEK293")
  }
  dge <- dge[, dge$samples$cell_line == line]

  # design matrix
  design <- model.matrix(~ treatment, data = dge$samples)
  colnames(design) <- gsub("^treatment", "", colnames(design))

  # filtering
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, ]

  # normalization
  dge$genes$spikein <- grepl("^ERCC-", row.names(dge))
  if (normalize == "global") {
    dge_global <- edgeR::calcNormFactors(dge[!dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_global$samples$norm.factors
  } else {
    dge_spikein <- edgeR::calcNormFactors(dge[dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_spikein$samples$norm.factors
  }

  # voom transformation
  v <- limma::voom(
    counts = dge,
    design = design,
    plot = FALSE,
    save.plot = FALSE
  )
  if (output == "v") {
    return(v)
  }

  # fit linear model
  fit <- limma::lmFit(object = v[, row.names(design)], design = design)
  fit <- limma::eBayes(fit)
  return(fit)
}


#' Fit a mixed linear model to RNA-seq data from fibroblasts
#'
#' Differential expression analysis using functions and methods from the
#' `edgeR` and `limma` packages.
#'
#' 1. The [edgeR::filterByExpr()] method is called to determine which genes have
#'    sufficiently large counts to be retained in a statistal analysis, based on
#'    the experimental design.
#' 2. The [edgeR::calcNormFactors()] method is used to calculate a TMM scaling
#'    factor, either based on all endogenous genes or only on the external ERCC
#'    spike-ins retained after filtering, as determined by the `normalize`
#'    argument.
#' 3. The [limma::voom()] function is used to transform raw count data to
#'    log2-counts per million (logCPM), estimate the mean-variance
#'    relationship and compute appropriate observation-level weights.
#' 4. The [limma::duplicateCorrelation()] function is used to estimate the
#'    correlation between between technical replicates (e.g. within each
#'    cell line).
#' 5. The [limma::voom()] function is called a second time, this time also
#'    providing the average estimated inter-duplicate correlation determined
#'    above.
#' 6. A linear model is fit for each gene with the [limma::lmFit()] function.
#' 7. Moderated t-statistics and log-odds of differential expression are
#'    calculated with the [limma::eBayes()] function using an empirical Bayes
#'    approach.
#' A linear mixed model is fit by REML individually for each gene, taking
#' account the correlation between technical replicates performed within each
#' cell line.
#'
#' @inheritParams fit_hek293
#' @importFrom limma voom duplicateCorrelation lmFit eBayes
#' @importFrom edgeR filterByExpr calcNormFactors [.DGEList
#' @return [limma::MArrayLM] fit object
#' @export
#' @examples
#' library(limma)
#' # global shift is detected with spike-in normalization
#' fit_spikein <- fit_fibroblasts(normalize = "spike-in")
#' limma::plotMA(fit_spikein, coef = "SUPT4H1",
#'     main = "Fibroblasts: SUPT4H1 RNAi vs Scramble (spike-in normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
#' # but not with global normalization
#' fit_global <- fit_fibroblasts(normalize = "global")
#' limma::plotMA(fit_global, coef = "SUPT4H1",
#'     main = "Fibroblasts: SUPT4H1 RNAi vs Scramble (global normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
fit_fibroblasts <- function(dge, normalize = c("spike-in", "global"),
                            output = c("fit", "v")) {
  normalize <- match.arg(normalize)
  output <- match.arg(output)

  if (missing(dge)) {
    dge <- load_dataset(experiment = "Fibroblasts")
  }

  # design matrix
  design <- model.matrix(~ rnai_treatment, data = dge$samples)
  colnames(design) <- gsub("^rnai_treatment", "", colnames(design))

  # filtering
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, ]

  # normalization
  dge$genes$spikein <- grepl("^ERCC-", row.names(dge))
  if (normalize == "global") {
    dge_global <- edgeR::calcNormFactors(dge[!dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_global$samples$norm.factors
  } else {
    dge_spikein <- edgeR::calcNormFactors(dge[dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_spikein$samples$norm.factors
  }

  # voom transformation (without correlation information)
  v <- voom(
    counts = dge,
    design = design,
    plot = FALSE,
    save.plot = FALSE
  )

  # Model duplicate correlation on previously voom-transformed data
  blocks <- with(v$targets, paste(rnai_treatment, replicate))
  corfit <- suppressWarnings(
    duplicateCorrelation(v[v$genes$spikein == FALSE, ],
                         design, block = blocks)
  )

  # secound round of voom transformation, now with correlation information
  v <- voom(
    counts = dge,
    design = design,
    correlation = corfit$consensus,
    plot = FALSE,
    save.plot = FALSE
  )
  if (output == "v") {
    return(v)
  }
  # fit linear model
  fit <- lmFit(object = v[, row.names(design)], design = design,
               block = blocks, correlation = corfit$consensus)
  fit <- eBayes(fit)
  return(fit)
}

#' Modeling the SUPT4H1 in fibroblasts effect separately for each siRNA conc.
#'
#' Differential expression analysis using functions and methods from the
#' `edgeR` and `limma` packages':
#' 1. The [edgeR::filterByExpr()] method is called to determine which genes have
#'    sufficiently large counts to be retained in a statistal analysis, based on
#'    the experimental design.
#' 2. The [edgeR::calcNormFactors()] method is used to calculate a TMM scaling
#'    factor, either based on all endogenous genes or only on the external ERCC
#'    spike-ins retained after filtering, as determined by the `normalize`
#'    argument.
#' 3. The [limma::voom()] function is used to transform raw count data to
#'    log2-counts per million (logCPM), estimate the mean-variance
#'    relationship and compute appropriate observation-level weights.
#' 4. A linear model is fit for each gene with the [limma::lmFit()] function.
#' 5. Moderated t-statistics and log-odds of differential expression are
#'    calculated with the [limma::eBayes()] function using an empirical Bayes
#'    approach.
#' A linear mixed model is fit by REML individually for each gene, taking
#' account the correlation between technical replicates performed within each
#' cell line.
#'
#' @inheritParams fit_hek293
#' @importFrom limma voom duplicateCorrelation lmFit eBayes topTable
#' @importFrom edgeR filterByExpr calcNormFactors [.DGEList
#' @importFrom purrr map map_df
#' @importFrom dplyr recode mutate
#' @importFrom tidyr replace_na
#' @importFrom rlang !!!
#' @export
#' @return A list of [limma::MArrayLM] fit objects, one for each siRNA
#' concentration.
fit_fibroblast_conc <- function(dge, normalize = c("spike-in", "global")) {
  normalize <- match.arg(normalize)

  if (missing(dge)) {
    dge <- load_dataset(experiment = "Fibroblasts")
  }

  fits <- purrr::map(
    setNames(c(5, 10),
             c("SUPT4H1_5nM", "SUPT4H1_10nM")),
    function(concentration) {
      dge_subset <- dge[, dge$samples$rnai_treatment == "scramble" |
                          dge$samples$treatment_conc == concentration]
      dge_subset$samples <- droplevels(dge_subset$samples)

      # design matrix
      design <- model.matrix(~ rnai_treatment, data = dge_subset$samples)
      colnames(design) <- gsub("^rnai_treatment", "", colnames(design))

      # filtering
      keep <- filterByExpr(dge_subset, design)
      dge_subset <- dge_subset[keep, ]

      # spike-in normalization
      dge_spikes <- calcNormFactors(dge_subset[dge_subset$genes$spikein, ])
      dge_subset$samples$norm.factors <- dge_spikes$samples$norm.factors

      # voom transformation
      v <- voom(
        counts = dge_subset,
        design = design,
        plot = FALSE,
        save.plot = FALSE
      )

      fit <- lmFit(object = v[, row.names(design)], design = design)
      fit <- eBayes(fit)
      return(fit)
    })
  return(fits)
}

#' Fit a mixed linear model to RNA-seq data from HEK293_71D cells
#'
#' Differential expression analysis using functions and methods from the
#' `edgeR` and `limma` packages.
#'
#' 1. The [edgeR::filterByExpr()] method is called to determine which genes have
#'    sufficiently large counts to be retained in a statistal analysis, based on
#'    the experimental design.
#' 2. The [edgeR::calcNormFactors()] method is used to calculate a TMM scaling
#'    factor, either based on all endogenous genes or only on the external ERCC
#'    spike-ins retained after filtering, as determined by the `normalize`
#'    argument.
#' 3. The [limma::voom()] function is used to transform raw count data to
#'    log2-counts per million (logCPM), estimate the mean-variance
#'    relationship and compute appropriate observation-level weights.
#' 4. A linear model is fit for each gene with the [limma::lmFit()] function.
#' 5. Moderated t-statistics and log-odds of differential expression are
#'    calculated with the [limma::eBayes()] function using an empirical Bayes
#'    approach.
#'
#' @param dge DGEList with raw counts
#' @param normalize Scalar string, calculate a TMM scaling factor using either
#' `spike-in` or `global` normalization.
#' @param output Scalar character specifying which output to return, either `v`
#' to return output of [limma::voom()] or `fit` to return the output of
#' [limma::eBayes()].
#' @importFrom limma voom lmFit eBayes
#' @importFrom edgeR filterByExpr calcNormFactors [.DGEList
#' @return [limma::MArrayLM] fit object
#' @export
#' @examples
#' library(limma)
#' # global shift is detected with spike-in normalization
#' fit_spikein <- fit_hek293_71D(normalize = "spike-in")
#'
#' limma::plotMA(fit_spikein, coef = "SUPT4H1",
#'     main = "HEK293 7D: SUPT4H1 RNAi vs Scramble (spike-in normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
#' # but not with global normalization
#' fit_global <- fit_hek293_71D(normalize = "global")
#' limma::plotMA(fit_global, coef = "SUPT4H1",
#'     main = "HEK293: SUPT4H1 RNAi vs Scramble (global normalization)",
#'     bg.col = adjustcolor("grey", 0.3), bg.cex = 0.7,
#'     status = factor(fit_spikein$genes$ercc_logFC),
#'     value = c("-2", "0", "0.58", "1"))
#' abline(h = 0, lty = 2)
fit_hek293_71D <- function(dge, normalize = c("spike-in", "global"),
                       output = c("fit", "v")) {
  normalize <- match.arg(normalize)
  output <- match.arg(output)

  if (missing(dge)) {
    dge <- load_dataset(experiment = "HEK293_71D")
  }

  # design matrix
  design <- model.matrix(~ rnai_treatment, data = dge$samples)
  colnames(design) <- gsub("^rnai_treatment", "", colnames(design))

  # filtering
  keep <- edgeR::filterByExpr(dge, design)
  dge <- dge[keep, ]

  # normalization
  dge$genes$spikein <- grepl("^ERCC-", row.names(dge))
  if (normalize == "global") {
    dge_global <- edgeR::calcNormFactors(dge[!dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_global$samples$norm.factors
  } else {
    dge_spikein <- edgeR::calcNormFactors(dge[dge$genes$spikein, ])
    dge$samples$norm.factors <- dge_spikein$samples$norm.factors
  }

  # voom transformation (without correlation information)
  v <- limma::voom(
    counts = dge,
    design = design,
    plot = FALSE,
    save.plot = FALSE
  )
  if (output == "v") {
    return(v)
  }

  # fit linear model
  fit <- limma::lmFit(object = v[, row.names(design)], design = design)
  fit <- limma::eBayes(fit)
  return(fit)
}
