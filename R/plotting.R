#' Generate an MA plot from a topTable
#'
#' @param tt, data.frame, typically the output of the [limma::topTable()]
#' function
#' @param xlim Numeric vector of length two that specifies the x-axis limits
#' @param ylim Numeric vector of length two that specifies the y-axis limits
#' @param nudge_spt4 Scalar numeric, the y-axis offset of the SPT4H1 label
#' @return A ggplot object
#' @importFrom ggplot2 ggplot scale_fill_continuous scale_color_brewer xlab ylab
#' theme element_blank aes stat_density2d geom_point geom_text geom_hline
#' coord_cartesian  scale_y_continuous scale_x_continuous guides guide_legend
#' @importFrom dplyr filter mutate %>%
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_ma <- function(tt, xlim = c(-3, 12), ylim = c(-6, 4), nudge_spt4 = 0.5) {
  ercc_medians <- calculate_ercc_medians(tt)
  ggplot(tt, aes(x = AveExpr, y = logFC)) +
    scale_fill_continuous(low = "white", high = "dodgerblue4", guide = FALSE) +
    scale_color_brewer(palette = "Set1", name = "ERCC logFC") +
    xlab("Average expression (log2)") +
    ylab(bquote(
      atop(
        italic("SUPT4H1")~"vs Scrambled RNAi",
        "log2 fold change")
    )) +
    theme(panel.grid = element_blank()) +
    stat_density2d(aes(fill = ..density..^0.25), geom = "tile",
                   contour = FALSE, n = 200) +
    geom_point(
      data = dplyr::filter(tt, hgnc_symbol == "SUPT4H1"),
      color = "black",
      size = 3
    ) +
    geom_text(
      data = dplyr::filter(tt, hgnc_symbol == "SUPT4H1") %>%
        dplyr::mutate(hgnc_symbol = sprintf("italic('%s')", hgnc_symbol)),
      aes(label = hgnc_symbol), color = "black", size  = 4,
      nudge_y = nudge_spt4, parse = TRUE) +
    geom_point(data = tt %>% filter(spikein),
               aes(color = factor(ercc_logFC,
                                  levels = c("-2", "0", "0.58", "1")))) +
    geom_hline(data = ercc_medians,
               aes(yintercept = median_fc,
                   color = factor(ercc_logFC,
                                  levels = c("-2", "0", "0.58", "1")))) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    guides(color = guide_legend(override.aes = list(fill = NA,
                                                    linetype = 0))) +
    theme(legend.position = "right")
}

#' Generate a boxplot from a topTable
#'
#' @param tt, data.frame, typically the output of the [limma::topTable()]
#' function
#' @param colors Character vector with colors for each box
#' @param ylim Numeric vector of length two that specifies the y-axis limits
#' @param nudge_x Scalar numeric, the x-axis offset for the labels
#' @param nudge_y Scalar numeric, the y-axis offset for the labels
#' @param nudge_spt4 Scalar numeric, the y-axis offset of the SPT4H1 label
#' @param cex Scalar numeric, character expansion for the SUPT4H1 shape
#' @param show_spt4 Scalar logical, highlight the log-fold change for SUPT4H1?
#' @return A ggplot object
#' @importFrom ggplot2 ggplot aes geom_hline geom_boxplot scale_y_continuous
#' scale_fill_manual scale_color_manual coord_cartesian xlab ylab facet_grid
#' theme element_blank geom_text geom_point
#' @importFrom dplyr mutate filter group_by summarise %>%
#' @importFrom tidyr replace_na
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_boxes <- function(tt, colors = brewer.pal(n = 8, "Set1"),
                       ylim = c(-3, 2), labels = FALSE,
                       nudge_x = 0.2,
                       nudge_y = 0.25,
                       nudge_spt4 = 0.25,
                       cex = 4,
                       show_spt4 = FALSE) {
  tt <- tt %>%
    dplyr::mutate(
      class = ifelse(spikein, "Spike-in controls", "Cellular\nRNA"),
      ercc_logFC = replace_na(ercc_logFC, "SUPT4H1 RNAi"),
      ercc_logFC = factor(ercc_logFC,
                          levels = c("-2", "0", "0.58", "1", "SUPT4H1 RNAi")))
  spt4 <- tt %>%
    dplyr::filter(hgnc_symbol == "SUPT4H1") %>%
    dplyr::mutate(hgnc_symbol = sprintf("italic('%s')", hgnc_symbol))
  median_logFC <- tt %>%
    dplyr::group_by(class, ercc_logFC) %>%
    dplyr::summarise(median_logFC = median(logFC))
  p <- tt %>%
    ggplot(aes(x = ercc_logFC,
               y = logFC,
               color = ercc_logFC,
               fill = ercc_logFC)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_boxplot(outlier.color = NA, color = "black") +
    scale_y_continuous(breaks = seq(-10, 10, by = 1)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    coord_cartesian(ylim = ylim) +
    xlab("Expected fold-change (log2)") +
    ylab("Observed fold change (log2)") +
    facet_grid(~ class, space = "free_x", scales = "free_x", drop = TRUE) +
    theme(panel.grid = element_blank(), legend.position = "none")
  if (labels == TRUE) {
    p <- p + geom_text(data = median_logFC,
                       aes(x = ercc_logFC,
                           y = median_logFC,
                           color = ercc_logFC,
                           label = round(median_logFC, 1)
                       ), nudge_y = nudge_y, nudge_x = nudge_x, size = cex)
  }
  if (show_spt4 == TRUE) {
    p <- p + geom_point(data = spt4, color = "black", size = 3) +
      geom_text(
        data = spt4,
        aes(label = hgnc_symbol), color = "black", size  = 4, parse = TRUE,
        nudge_y = nudge_spt4)

  }
  return(p)
}

#' Plot the normalized CPMs for a single gene symbol
#'
#' @param gene_symbol An HGNC gene symbol, must match an entry in
#' `dge$genes$hgnc_symbol`.
#' @param dge A DGEList list object with appropriate normalization applied
#' @param labels A named character vector that will be used to rename elements
#' of the `dge$samples$treatment` and `dge$samples$cell_line` columns.
#' @return A ggplot object
#' @importFrom edgeR cpm
#' @importFrom ggplot2 ggplot aes geom_jitter facet_grid scale_fill_manual
#' scale_shape theme ylab xlab ggtitle theme element_blank
#' @importFrom dplyr mutate left_join recode %>%
#' @importFrom tidyr gather
#' @importFrom stringr str_split_fixed
#' @importFrom rlang !!!
#' @export
plot_expression <- function(gene_symbol, dge, labels) {
  stopifnot(gene_symbol %in% dge$genes$hgnc_symbol)
  edgeR::cpm(dge[gene_symbol, ], normalized.lib.sizes = TRUE,
             log = FALSE) %>%
    as.data.frame(stringsAsFactor = FALSE) %>%
    tidyr::gather(key = "sample_id", value = "cpm") %>%
    dplyr::mutate(hgnc_symbol = gene_symbol) %>%
    dplyr::left_join(
      dplyr::mutate(
        dge$samples,
        replicate = paste("Replicate",
                          str_split_fixed(name, pattern = "_", n = 3)[, 3],
                          sep = "_"),
        treatment = recode(treatment, !!!labels),
        cell_line = recode(cell_line, !!!labels)
      ),
      by = "sample_id") %>%
    ggplot(aes(x = treatment, y = cpm, fill = treatment)) +
    geom_jitter(width = 0.1, height = 0, size = 4, alpha = 0.7, shape = 21,
                color = "grey") +
    facet_grid(. ~ cell_line, scales = "free_y") +
    scale_fill_manual(values = c(
      SUPT4H1 = "orange", Scramble = "darkgrey",
      "Scramble\nsiRNA" = "darkgrey"),
      name = "", guide = FALSE) +
    scale_shape(name = "") +
    theme(legend.position = "bottom") +
    ylab("Normalized CPM") +
    xlab("") +
    ggtitle(gene_symbol) +
    theme(panel.grid = element_blank())
}

#' Plot SUPT4H1 RNAi effect by biotype
#'
#' @param tt, data.frame, typically the output of the [limma::topTable()]
#' function
#' @importFrom dplyr filter group_by mutate ungroup %>%
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_manual geom_hline xlab
#' ylab coord_flip theme element_blank
#' @return ggplot object
#' @export
plot_biotypes <- function(tt) {
  tt %>% dplyr::filter(
    !gene_biotype %in% c("control", "rRNA")) %>%
    dplyr::group_by(gene_biotype) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      protein_coding = gene_biotype == "protein_coding",
      biotype = sprintf("%s (%s)", gsub("_", " ", gene_biotype), n),
      biotype = factor(
        biotype, levels = rev(sort(unique(biotype))))) %>%
    ggplot(aes(x = biotype, y = logFC, fill = protein_coding)) +
    geom_boxplot(outlier.colour = NA) +
    scale_fill_manual(values = c("TRUE" = "orange", "FALSE" = "grey"),
                      name = "", guide = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("") +
    ylab(bquote(
      atop(
        italic("SUPT4H1")~"vs Scrambled RNAi",
        "log2 fold change")
    )) +
    coord_flip(ylim = c(-4, 2)) +
    theme(panel.grid = element_blank())
}

#' Create density plots of flow cytometry results
#'
#' @param x data.frame with two columns:
#' - signal
#' - condition
#' @param red Scalar character, the condition displayed in red.
#' @param grey Scalar character, the condition displayed in grey.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_density ylab geom_point scale_color_manual
#' scale_x_continuous guides guide_legend coord_cartesian theme
#' @export
#' @examples
#' plot_flow_density(load_flow_data("siRNA"))
plot_flow_density <- function(x, red = "SUPT4H1", grey="Scramble") {
  x[x$signal < 0, "signal"] <- 0L
  x_red <- dplyr::filter(x, condition == red)
  x_grey <- dplyr::filter(x, condition == grey)
  ggplot(x, aes(x = signal + 1)) +
    geom_density(data = x_grey, color = "darkgrey", fill = "grey", size = 0.7) +
    geom_density(data = x_red, color = "firebrick", size = 0.9) +
    ylab("Density") +
    theme(panel.grid = element_blank()) +
    # add invisible points to create a custom legend
    geom_point(
      data = data.frame(class = c(red, grey),
                        x = 0, y = 0),
      aes(x = x, y = y, color = class), shape = NA) +
    scale_color_manual(
      values = setNames(c("darkgrey", "firebrick"), c(grey, red)),
      name = "") +
    scale_x_continuous(trans = scales::log1p_trans(),
                       breaks = 10^seq(1, 4, by = 1)) +
    guides(color = guide_legend(
      override.aes = list(fill = NA, shape = 15, linetype = 0, size = 3))) +
    coord_cartesian(xlim = c(10, 1e4)) +
    xlab("Fluorescence (FITC, A.U.)") +
    theme(legend.position = "bottom")
}

#' Plot differential expression statistics for a gene set
#'
#' @param tt, data.frame, typically the output of the [limma::topTable()]
#' @param gene_sets a list of one or more gene sets, each with ensembl
#' identifiers
#' @param xlim Numeric vector of length two, the x-axis limits
#' @param y Scalar character specifying which statistic to plot.
#' @return ggplot object
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_density geom_rug geom_vline xlab ylab
#' coord_cartesian theme element_blank facet_wrap geom_point scale_color_manual
#' guides guide_legend aes_string
plot_gsea <- function(tt, gene_sets, id, xlim = c(-3, 3), y = c("logFC", "y")) {
  y <- match.arg(y)
  x_lab <- switch(
    y, "logFC" = bquote(
      atop(
        italic("SUPT4H1")~"vs Scrambled RNAi",
        "log2 fold change")
    ),
    "t" = bquote(
      atop(
        italic("SUPT4H1")~"vs Scrambled RNAi",
        "t-statistic")
    )
  )
  genes <- gene_sets[[id]]
  tt_react <- dplyr::filter(tt, ensembl_gene_id %in% genes)
  main <- paste(strwrap(paste0(toupper(substr(id, 1, 1)),
                               substr(id, 2, nchar(id))),
                        40), collapse = "\n")
  p <- tt %>%
    ggplot(aes_string(x = y)) +
    geom_density(color = "darkgrey", fill = "grey", size = 0.7) +
    geom_density(data = tt_react, color = "firebrick", size = 0.9) +
    geom_rug(data = tt_react, color = "firebrick") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    coord_cartesian(xlim = xlim) +
    xlab(x_lab) +
    ylab("Density") +
    theme(panel.grid = element_blank()) +
    ggtitle(main)
  if ("group" %in% colnames(tt)) {
    p <- p + facet_wrap(~ group)
  }
  # add invisible points to create a custom legend
  p <- p + geom_point(
    data = data.frame(class = c("All genes", "Gene set"),
                      x = 0, y = 0),
    aes(x = x, y = y, color = class), shape = NA) +
    scale_color_manual(
      values = c("All genes" = "darkgrey", "Gene set" = "firebrick"),
      name = "") +
    guides(color = guide_legend(
      override.aes = list(fill = NA,shape = 15, linetype = 0, size = 3)))
  return(p)
}

#' Plot normalized CPMs for a single human gene
#'
#' @param dge A DGEList list object with appropriate normalization applied
#' @param gene_symbol Scalar character, a human gene identifier
#' @return ggplot object
#' @export
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom rlang ensym
#' @importFrom ggplot2 ggplot xlab ylab ggtitle geom_jitter scale_fill_manual
#' scale_shape theme theme_linedraw
#' @importFrom edgeR cpm
plot_scatter <- function(dge, gene_symbol = "GFP", x = "rnai_treatment") {
  x <- rlang::ensym(x)
  edgeR::cpm(dge[dge$genes$hgnc_symbol == gene_symbol, ],
             normalized.lib.sizes = TRUE,
             log = FALSE) %>%
    as.data.frame(stringsAsFactor = FALSE) %>%
    tidyr::gather(key = "sample_id", value = "cpm") %>%
    dplyr::left_join(
        dge$samples %>%
          as.data.frame() %>%
          tibble::rownames_to_column("sample_id"),
        by = "sample_id") %>%
    ggplot(aes(x = !!x, y = cpm, fill = !!x)) +
    geom_jitter(width = 0.1, height = 0, size = 4, alpha = 0.7, shape = 21,
                color = "grey") +
    scale_fill_manual(values = c(
      SUPT4H1 = "orange", Scramble = "darkgrey", scramble = "darkgrey"),
      name = "", guide = FALSE) +
    scale_shape(name = "") +
    ylab("Normalized CPM") +
    xlab("") +
    ggtitle(gene_symbol) +
    theme_linedraw() +
    theme(panel.grid = element_blank(), legend.position = "bottom")
}

