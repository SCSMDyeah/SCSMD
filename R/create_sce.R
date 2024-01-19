#' Create SingleCellExperiment object from count data
#'
#' @param counts gene-by-cell count matrix
#' @param col.data a dataframe containing cell informations
#' @param row.data a dataframe containing gene informations
#' @param scale.factor scalar used for per cell sum of count normalization
#'
#' @return a SingleCellExperiment object
#' @importFrom utils installed.packages
#'
create_sce_from_counts <- function(counts, col.data, row.data = NULL, scale.factor = 10000) {

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
    sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts)),
                                                         colData = col.data)
  } else {
    sceset <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts)),
                                                         colData = col.data, rowData = row.data)
  }
  # this function writes to logcounts slot exprs(sceset) <-
  # log2(calculateCPM(sceset, use_size_factors = FALSE) + 1)
  SingleCellExperiment::logcounts(sceset) = log2(t(t(counts)/colSums(counts)) *
                                                   scale.factor + 1)
  # use gene names as feature symbols
  SummarizedExperiment::rowData(sceset)$feature_symbol = rownames(sceset)
  # remove features with duplicated names
  if (is.null(row.data)) {
    sceset <- sceset[!duplicated(SummarizedExperiment::rowData(sceset)$feature_symbol),
    ]
  }
  # QC isSpike(sceset, 'ERCC') <- grepl('^ERCC-', rownames(sceset)) sceset <-
  # calculateQCMetrics(sceset, feature_controls = list('ERCC' = isSpike(sceset,
  # 'ERCC')))
  return(sceset)
}

create_sce_from_logcounts <- function(logcounts, colData, rowData = NULL) {
  if(is.null(rowData)) {
    sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)),
                                   colData = colData)
  } else {
    sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)),
                                   colData = colData,
                                   rowData = rowData)
  }
  # use gene names as feature symbols
  rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if(is.null(rowData)) {
    sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
  }
  # QC
  isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
  return(sceset)
}
