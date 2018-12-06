#' Concentrations for ERCC spike-in mixes
#'
#' @export
#' @return data.frame with one row for each ERCC spike-in transcript
#' @examples
#' head(ercc())
ercc <- function() {
  read.csv(system.file("extdata", "ERCC.csv", package = "SUPT4H1"),
           stringsAsFactors = FALSE)
}

#' Ensembl identifiers for Peptide Chain Elongation Reactome gene set
#'
#' Entrez gene identifiers from the `Peptide Chain Elongation Reactome` gene
#' set were curated from the `reactome.db` Bioconductor package
#' (version 1.64.0) and mapped to ensembl gene identifiers using the
#' `org.Hs.eg.db` Bioconductor package (version 3.6.0).
#'
#' @return list with a single element, a named vector with ensembl identifiers
#' annotated in the `Peptide Chain Elongation` Reactome gene set.
#' @export
#' @examples
#' peptide_chain_elongation()
peptide_chain_elongation <- function() {
  list(
    "Peptide Chain Elongation" = setNames(
      c(
      "ENSG00000197728", "ENSG00000136942", "ENSG00000163923",
      "ENSG00000280969", "ENSG00000165496", "ENSG00000156508",
      "ENSG00000167658", "ENSG00000167658", "ENSG00000163584",
      "ENSG00000142541", "ENSG00000130255", "ENSG00000168028",
      "ENSG00000198755", "ENSG00000185088", "ENSG00000037241",
      "ENSG00000100316", "ENSG00000140986", "ENSG00000174444",
      "ENSG00000122406", "ENSG00000089009", "ENSG00000147604",
      "ENSG00000148303", "ENSG00000161016", "ENSG00000163682",
      "ENSG00000147403", "ENSG00000142676", "ENSG00000197958",
      "ENSG00000167526", "ENSG00000174748", "ENSG00000265681",
      "ENSG00000063177", "ENSG00000105640", "ENSG00000108298",
      "ENSG00000122026", "ENSG00000116251", "ENSG00000198242",
      "ENSG00000114391", "ENSG00000161970", "ENSG00000131469",
      "ENSG00000156482", "ENSG00000166441", "ENSG00000108107",
      "ENSG00000162244", "ENSG00000071082", "ENSG00000144713",
      "ENSG00000109475", "ENSG00000182899", "ENSG00000165502",
      "ENSG00000145592", "ENSG00000197756", "ENSG00000172809",
      "ENSG00000198918", "ENSG00000229117", "ENSG00000241343",
      "ENSG00000089157", "ENSG00000137818", "ENSG00000177600",
      "ENSG00000140988", "ENSG00000149273", "ENSG00000145425",
      "ENSG00000198034", "ENSG00000129824", "ENSG00000083845",
      "ENSG00000137154", "ENSG00000171863", "ENSG00000142937",
      "ENSG00000170889", "ENSG00000124614", "ENSG00000142534",
      "ENSG00000112306", "ENSG00000110700", "ENSG00000164587",
      "ENSG00000115268", "ENSG00000134419", "ENSG00000105193",
      "ENSG00000182774", "ENSG00000231500", "ENSG00000105372",
      "ENSG00000008988", "ENSG00000171858", "ENSG00000186468",
      "ENSG00000138326", "ENSG00000118181", "ENSG00000197728",
      "ENSG00000177954", "ENSG00000143947", "ENSG00000233927",
      "ENSG00000213741", "ENSG00000221983", "ENSG00000188846",
      "ENSG00000125691"),
      c("RPS26", "RPL35", "RPL39L", "RPS4Y2", "RPL10L", "EEF1A1", "EEF2",
        "EEF2", "RPL22L1", "RPL13A", "RPL36", "RPSA", "RPL10A", "RPS27",
        "RPL26L1", "RPL3", "RPL3L", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A",
        "RPL8", "RPL9", "RPL10", "RPL11", "RPL12", "RPL13", "RPL15", "RPL17",
        "RPL18", "RPL18A", "RPL19", "RPL21", "RPL22", "RPL23A", "RPL24",
        "RPL26", "RPL27", "RPL30", "RPL27A", "RPL28", "RPL29", "RPL31", "RPL32",
        "RPL34", "RPL35A", "RPL36AL", "RPL37", "RPL37A", "RPL38", "RPL39",
        "RPL41", "RPL36A", "RPLP0", "RPLP1", "RPLP2", "RPS2", "RPS3", "RPS3A",
        "RPS4X", "RPS4Y1", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10",
        "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17",
        "RPS18", "RPS19", "RPS20", "RPS21", "RPS23", "RPS24", "RPS25", "RPS26",
        "RPS27", "RPS27A", "RPS28", "RPS29", "UBA52", "RPL14", "RPL23")
    )
  )
}
