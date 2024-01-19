#' Create SingleCellExperiment object from normalized data
#'
#' @param normcounts gene-by-cell matrix of normalized expressions
#' @param col.data a dataframe containing cell informations
#' @param row.data a dataframe containing gene informations
#'
#' @return a SingleCellExperiment object
#' @importFrom SingleCellExperiment SingleCellExperiment rowData
#' @importFrom utils installed.packages
create_sce_from_normcounts <- function(normcounts, col.data, row.data = NULL) {

  if (!"SingleCellExperiment" %in% utils::installed.packages()) {
    BiocManager::install("SingleCellExperiment")
    if (!"SingleCellExperiment" %in% utils::installed.packages()) {
      stop("Please install SingleCellExperiment r package!")
    }
  }

  if (!"SummarizedExperiment" %in% utils::installed.packages()) {
    BiocManager::install("SummarizedExperiment")
    if (!"SummarizedExperiment" %in% utils::installed.packages()) {
      stop("Please install SummarizedExperiment r package!")
    }
  }

  if (is.null(row.data)) {
    sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)),
                                                         colData = col.data)
  } else {
    sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)),
                                                         colData = col.data, rowData = row.data)
  }
  SingleCellExperiment::logcounts(sceset) <- SingleCellExperiment::normcounts(sceset)
  # use gene names as feature symbols
  SummarizedExperiment::rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if (is.null(rowData)) {
    sceset <- sceset[!duplicated(SummarizedExperiment::rowData(sceset)$feature_symbol),
    ]
  }
  # QC isSpike(sceset, 'ERCC') <- grepl('^ERCC-', rownames(sceset))
  return(sceset)
}
