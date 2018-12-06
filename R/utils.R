#' Format numbers in a human-friendly format as axis labels
#'
#' @param x Vector of numbers to be reformatted
#' @param signif Scalar integer, the number of significant digits to include
#' @export
#' @return Character vector
#' @examples
#' human_numbers(c(-100, 1, 1000, 1.5e6))
human_numbers <- function(x, signif = 1){
  humanity <- function(y){

    if (!is.na(y)) {
      tn <- round(abs(y) / 1e12, signif)
      b <- round(abs(y) / 1e9, signif)
      m <- round(abs(y) / 1e6, signif)
      k <- round(abs(y) / 1e3, signif)

      if (y >= 0) {
        y_sign <- ""
      } else {
        y_sign <- "-"
      }

      if ( k < 1 ) {
        paste0(y_sign, round(abs(y), signif ))
      } else if (m < 1) {
        paste0(y_sign,  k , "k")
      } else if (b < 1) {
        paste0(y_sign, m ,"m")
      }else if (tn < 1) {
        paste0(y_sign, b ,"bn")
      } else {
        paste0(y_sign,  comma(tn), "tn")
      }
    } else if (is.na(y) | is.null(y)) {
      "-"
    }
  }
  vapply(x, humanity, character(1))
}

#' Calculate the median expression for each ERCC group
#'
#' @param tt, data.frame, typically the output of the [limma::topTable()]
#' function
#' @return A tibble with one row for each ERCC spike-in group
#' @importFrom dplyr group_by summarize filter %>%
#' @export
calculate_ercc_medians <- function(tt) {
  if ("cell_line" %in% colnames(tt)) {
    ercc_medians <- tt %>%
      dplyr::group_by(cell_line, ercc_logFC) %>%
      dplyr::summarize(median_fc = median(logFC)) %>%
      dplyr::filter(!is.na(ercc_logFC))
  } else {
    ercc_medians <- tt %>%
      dplyr::group_by(ercc_logFC) %>%
      dplyr::summarize(median_fc = median(logFC)) %>%
      dplyr::filter(!is.na(ercc_logFC))
  }
  return(ercc_medians)
}
