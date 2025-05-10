#' Baron et al. (2016) Human Pancreas Single-Cell RNA-Seq Reference Data
#'
#' A SingleCellExperiment object containing human pancreas single-cell RNA-seq data
#' from Baron et al. (2016). This dataset serves as the reference dataset in the
#' scLearn package for cell type annotation.
#'
#' @format A SingleCellExperiment object with 20125 genes (rows) and 8569 cells (columns).
#'   The object contains:
#'   \describe{
#'     \item{counts}{Raw count matrix}
#'     \item{logcounts}{Log-normalized expression values}
#'     \item{colData}{Cell metadata including:
#'       \itemize{
#'         \item{cell_type1: Cell type annotations}
#'         \item{other metadata fields...}
#'       }
#'     }
#'   }
#'
#' @source \url{https://doi.org/10.1016/j.cels.2016.08.011}
#'
#' @references
#' Baron, M., Veres, A., Wolock, S. L., Faust, A. L., Gaujoux, R., Vetere, A., ... & Melton, D. A. (2016).
#' A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell population structure.
#' Cell systems, 3(4), 346-360.
#'
#' @examples
#' library(SingleCellExperiment)
#' data(RefCellData)
#' # Access expression matrix
#' assay(RefCellData, "logcounts")[1:5,1:5]
#' # Access cell metadata
#' head(colData(RefCellData))
"RefCellData"

#' Muraro et al. (2016) Human Pancreas Single-Cell RNA-Seq Query Data
#'
#' A SingleCellExperiment object containing human pancreas single-cell RNA-seq data
#' from Muraro et al. (2016). This dataset serves as the query dataset in the
#' scLearn package for cell type annotation.
#'
#' @format A SingleCellExperiment object with 19127 genes (rows) and 2122 cells (columns).
#'   The object contains:
#'   \describe{
#'     \item{counts}{Raw count matrix}
#'     \item{logcounts}{Log-normalized expression values}
#'     \item{colData}{Cell metadata including:
#'       \itemize{
#'         \item{cell_type1: Cell type annotations}
#'         \item{other metadata fields...}
#'       }
#'     }
#'   }
#'
#' @source \url{https://doi.org/10.1016/j.cmet.2016.08.020}
#'
#' @references
#' Muraro, M. J., Dharmadhikari, G., Grün, D., Groen, N., Dielen, T., Jansen, E., ... & van Oudenaarden, A. (2016).
#' A single-cell transcriptome atlas of the human pancreas.
#' Cell systems, 3(4), 385-394.
#'
#' @examples
#' library(SingleCellExperiment)
#' data(QueryCellData)
#' # Access expression matrix
#' assay(QueryCellData, "logcounts")[1:5,1:5]
#' # Access cell metadata
#' head(colData(QueryCellData))
"QueryCellData"
