#' Return a list of ensembl identifiers for each Reactome pathway
#'
#' @param min_size Scalar integer, the minimum size of each set
#' @param max_size Scalar integer, the maximum size of each set
#' @importFrom reactome.db reactomePATHID2NAME reactomeEXTID2PATHID
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr inner_join
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
#' @examples
#' gene_sets <- reactome_sets()
#' gene_sets[1:2]
reactome_sets <- function(min_size = 10, max_size = Inf) {
  stopifnot(max_size > 0 & max_size > min_size)
  pathway.names <- as.data.frame(reactomePATHID2NAME)
  pathway.names <- pathway.names[
    grepl("Homo sapiens: ", pathway.names$path_name), ]
  pathway.names$path_name <- sub("Homo sapiens: ", "",
                                 pathway.names$path_name, fixed = TRUE)

  gene.pathway <- as.data.frame(reactomeEXTID2PATHID)[, c(2, 1)]
  gene.pathway <- gene.pathway[gene.pathway$DB_ID %in% pathway.names$DB_ID, ]
  gene.pathway$gene_id <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = gene.pathway$gene_id,
           column = "ENSEMBL", keytype = "ENTREZID"))
  gene.pathway <- na.omit(gene.pathway)
  gene_sets <- dplyr::inner_join(pathway.names, gene.pathway, by = "DB_ID") %>%
    split(f = .$path_name, x = .$gene_id)
  gene_sets[lengths(gene_sets) > min_size & length(gene_sets) < max_size]
}

#' Convert camera output into a datatable
#'
#' @param results data.frame, the output of the [limma::camera()] function
#' @export
#' @return a datatable htmlwidget
#' @importFrom dplyr arrange
#' @importFrom DT datatable formatSignif
camera_table <- function(results) {
  results %>%
    dplyr::arrange(Direction, FDR) %>%
    DT::datatable(
      colnames = c("Gene Set", "NGenes", "Direction", "p-Value", "FDR"),
      rownames = FALSE,
      filter = "bottom",
      extensions = c('Buttons'),
      options = list(
        columnDefs = list(list(className = 'dt-center', targets = 0:4)),
        pageLength = 10,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )
    ) %>%
    DT::formatSignif(digits = 3, columns = c("PValue", "FDR"))
}

