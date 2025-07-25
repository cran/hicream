if (getRversion() >= "2.15.1") utils::globalVariables("computeClustering")

#' @title Performs Constrained 2D Agglomerative Clustering
#'
#' @description This function performs a connectivity constrained 2D
#' agglomerative clustering using \code{scikit-learn} function
#' \code{AgglomerativeClustering} and outputs an object of class \code{hclust}
#' that stores the hierarchy of merges and value of criterion at each merge. It
#' also outputs the optimal level of the hierarchy with respect to the elbow
#' heuristic.
#'
#' @param counts an object of class
#' \code{\link[InteractionSet:interactions]{InteractionSet}} obtained from the
#' function \code{\link{loadData}}.
#' @param nbClust integer. Number of clusters to obtain. Set to \code{NULL} by
#' default.
#'
#' @return An object of class \code{res2D} containing:
#' \item{tree}{an object of class \code{hclust}}
#' \item{nbClust}{the number of clusters corresponding either to the value
#' passed by the user or to the optimal level of clusters as provided by the
#' elbow heuristic}
#' \item{clustering}{obtained clustering}
#'
#' @author Élise Jorge \email{elise.jorge@inrae.fr}\cr
#' Sylvain Foissac \email{sylvain.foissac@inrae.fr}\cr
#' Pierre Neuvial \email{pierre.neuvial@math.univ-toulouse.fr}\cr
#' Nathalie Vialaneix \email{nathalie.vialaneix@inrae.fr}
#'
#' @export
#'
#' @import reticulate
#' @importFrom Matrix sparseMatrix
#' @examples
#' data("pighic")
#' \donttest{res2D <- AggloClust2D(pighic$data)
#' if (!is.null(res2D)) {# in case Python or modules are not available
#'   clusters <- res2D$clustering
#'   print(res2D)
#'   summary(res2D)
#'   plot(res2D)
#' }}
#'
AggloClust2D <- function(counts, nbClust = NULL) {
  call <- sys.call()
  # transform data into count format
  counts <- fromDGE2Counts(counts)
  # order counts
  counts <- counts[with(counts, order(bin1, bin2)), ]
  # build neighborhood relative to connectivity constraint
  neighbors <- buildNeighborsR(counts)
  nbInt <- length(counts$bin1)
  # fill sparse neighboring matrix
  matNeighbors <- neighborsToMat(neighbors, nbInt)
  
  # check python and module availability
  modules_avail <- reticulate::py_available(initialize = TRUE) &&
    reticulate::py_module_available("sklearn") &&
    reticulate::py_module_available("kneebow") &&
    reticulate::py_module_available("pandas") &&
    reticulate::py_module_available("numpy")
  if (!modules_avail) {
    message("Python or required Python modules are not available. ",
            "Please install them to use this function.")
    return(NULL)
  }
  
  sourceFile <- file.path(system.file(package = "hicream"), "python", 
                          "2Dclust.py")
  source_python(sourceFile)
  if (is.null(nbClust)) {
    res2D <- computeClustering(counts, matNeighbors)
  } else {
    nbClust <- as.integer(nbClust)
    res2D <- computeClustering(counts, matNeighbors, nbClust = nbClust)
  }
  merge <- res2D[[1]]
  height <- res2D[[2]]
  nbClust <- res2D[[3]]
  hc <- fromAgglomerative2hclust(merge, height, call)
  outList <- list(hc, nbClust, res2D[[4]])
  names(outList) <- c("tree", "nbClust", "clustering")
  # Assign class to the list
  class(outList) <- "res2D"
  return(outList)
}

##### Methods for res2D object ####
#' @exportS3Method
#' @param x a \code{res2D} object to print
#' @param ... not used
#' @rdname AggloClust2D
print.res2D <- function(x, ...) {
  # print type of test
  cat("Tree obtained from constrained 2D clustering.\n")

  # print tree info
  print(x$tree)
  cat("\n")

  # print optimal number of clusters
  cat("Optimal number of clusters:", x$nbClust, "\n\n")

  # clustering
  cat("Clustering:\n")
  if (length(x$clustering) < 10) {
    cat(x$clustering)
  } else {
    cat(head(x$clustering, 10), "...\n")
  }
  invisible(NULL)
}

#' @exportS3Method
#' @param object a \code{res2D} object to summarize
#' @param ... not used
#' @rdname AggloClust2D
summary.res2D <- function(object, ...) {
  # print type of test
  cat("Summary of 2D constrained clustering results.")
  cat("\n\n")

  summary(object$tree)
}

#' @exportS3Method
#' @param x a \code{res2D} object to plot
#' @rdname AggloClust2D
plot.res2D <- function(x, ...) {
  plot(x$tree, hang = -1)
}
