#' @title Perform diffHic Test
#'
#' @description This function performs diffHic differential analysis for every
#'  pixel of the matrix.
#'
#' @param matrices an object of class
#' \code{\link[InteractionSet:interactions]{InteractionSet}} obtained from the
#' function \code{\link{loadData}}. Conditions correspond to the columns (one
#' column per replicate).
#' @param cond a vector indicating the condition of each column of
#' \code{matrices}.
#' @param outFile path to export outputs of the function, set to \code{NULL} by
#' default, in which case results are not exported.
#' @param filterLc logical. Whether to filter out low counts or not. Set to
#' \code{FALSE} by default, in which case no filtering is performed. See
#' \code{\link[diffHic:filterTrended]{filterTrended}} for more details.
#' @param filterFlank logical. Whether to filter out on enriched pairs. Set to
#' \code{FALSE} by default, in which case no filtering is performed. See
#' \code{\link[diffHic:enrichedPairs]{enrichedPairs}} and
#' \code{\link[diffHic:filterPeaks]{filterPeaks}} for more details.
#' @param flank flank parameter used only if \code{filterFlank = TRUE}. Set to
#' \code{NULL} by default. See
#' \code{\link[diffHic:enrichedPairs]{enrichedPairs}} for more details.
#'
#' @return An object of class \code{resdiff} with the following entries:
#' \item{region1}{the first bin of the interaction}
#' \item{region2}{the second bin of the interaction}
#' \item{p.value}{the p-value of the \code{\link[diffHic]{diffHic}} test}
#' \item{p.adj}{the adjusted p-value of the \code{\link[diffHic]{diffHic}} test}
#' \item{logFC}{the log2-fold-change of the interaction}
#'
#' @author Élise Jorge \email{elise.jorge@inrae.fr}\cr
#' Sylvain Foissac \email{sylvain.foissac@inrae.fr}\cr
#' Pierre Neuvial \email{pierre.neuvial@math.univ-toulouse.fr}\cr
#' Nathalie Vialaneix \email{nathalie.vialaneix@inrae.fr}
#'
#' @export
#'
#' @importFrom diffHic filterTrended
#' @importFrom diffHic enrichedPairs
#' @importFrom diffHic filterPeaks
#' @importFrom stats relevel
#' @importFrom stats model.matrix
#' @importFrom stats p.adjust
#' @importFrom stats relevel
#' @importFrom csaw asDGEList
#' @importFrom edgeR estimateDisp
#' @importFrom edgeR glmQLFTest
#' @importFrom edgeR glmQLFit
#' @importFrom adjclust plotSim
#' @importFrom viridis scale_fill_viridis
#'
#' @examples
#' data("pighic")
#' resdiff <- performTest(pighic$data, pighic$conditions)
#' resdiff
#' summary(resdiff)
#' plot(resdiff)
#' plot(resdiff, whichPlot = "p.adj")
#' plot(resdiff, whichPlot = "logFC")
#'
performTest <- function(matrices, cond, outFile = NULL, filterLc = FALSE,
                        filterFlank = FALSE, flank = NULL) {
  # Remove low counts.
  if (filterLc) { # based on background estimation and trended filtering strategy
    trended <- filterTrended(matrices)
    trendKeep <- trended$abundances > trended$threshold
    # trendKeep <- trendKeep & countKeep
    matrices <- matrices[trendKeep, ]
    colnames(matrices) <- NULL
  }
  # Filter based on enriched pairs.
  if (as.logical(filterFlank)) {
    enriched <- enrichedPairs(matrices, flank = flank)
    enrichments <- filterPeaks(enriched, get.enrich = TRUE)
    peakKeep <- filterPeaks(matrices, enrichments)
    matrices <- enriched[peakKeep, ]
  }

  # create groups
  colnames(matrices) <- paste(cond)
  cond <- factor(cond)
  # perform test
  design <- model.matrix(~cond)
  curMat <- asDGEList(matrices)
  curMat <- estimateDisp(curMat, design)
  result <- glmQLFit(curMat, design, robust = TRUE)
  result <- glmQLFTest(result, coef = 2)
  adjpval <- p.adjust(result$table$PValue, method = "BH")
  outSetsStart <- rowData(matrices)[, 1]
  outSetsEnd <- rowData(matrices)[, 2]

  # Store results in a list
  outTable <- data.frame(
    "region1" = outSetsStart,
    "region2" = outSetsEnd,
    "p.value" = result$table$PValue,
    "p.adj" = adjpval,
    "logFC" = result$table$logFC
  )

  # order bins
  outTable <- outTable[with(outTable, order(region1, region2)), ]

  # if the user specified an output destination, the output list is exported
  if (!is.null(outFile)) {
    write.csv(outTable, file = outFile)
  }

  # Assign class to the list
  class(outTable) <- "resdiff"

  return(outTable)
}


##### Methods for resdiff object ####
#' @exportS3Method
#' @param x a \code{resdiff} object to print
#' @param object a \code{resdiff} object to print
#' @param ... not used
#' @rdname performTest
print.resdiff <- function(x, ...) {
  # print type of test
  cat("Testing for differential interactions using diffHic.")
  cat("\n\n")

  # build data.frame
  dfPrint <- data.frame(region1 = x$region1, region2 = x$region2)
  dfPrint$p.value <- x$p.value
  dfPrint$p.adj <- x$p.adj
  dfPrint$logFC <- x$logFC

  # print first 10 rows
  print(dfPrint, max = 10)
}

#' @exportS3Method
#' @rdname performTest
summary.resdiff <- function(object, ...) {
  # print type of test
  cat("Summary of diffHic results.")
  cat("\n\n")
  dfRes <- data.frame(
    p.value = object$p.value, p.adj = object$p.adj,
    logFC = object$logFC
  )
  # Summarize only p.value, p.adj and logFC.
  summary(dfRes)
}

#' @exportS3Method
#' @param x a \code{resdiff} object to plot
#' @param whichPlot a character string indicating which plot to display.
#' Possible values are \code{"p.value"}, \code{"p.adj"} and \code{"logFC"}. Set
#' to \code{"p.value"} by default.
#' @rdname performTest
plot.resdiff <- function(x, whichPlot = c("p.value", "p.adj", "logFC"), ...) {
  whichPlot <- match.arg(whichPlot)
  uniqueBins <- unique(c(x$region1, x$region2))
  minBins <- min(uniqueBins)
  maxBins <- max(uniqueBins)
  p <- maxBins - minBins + 1
  if (whichPlot == "p.value") {
    mat <- matrix(1, p, p)
    mat[cbind(x$region1 - minBins + 1, x$region2 - minBins + 1)] <- x$p.value
    mat[cbind(x$region2 - minBins + 1, x$region1 - minBins + 1)] <- x$p.value
    allBreaks <- seq(-0.694, 0.3, by = 0.095)
    labelsPval <- round((exp(allBreaks) - 0.5), digits = 3)
    suppressMessages({
      pp <- plotSim(mat) + scale_fill_viridis(
        name = "p-value", breaks = allBreaks,
        labels = labelsPval, direction = -1
      )
    })
  }

  if (whichPlot == "p.adj") {
    mat <- matrix(1, p, p)
    mat[cbind(x$region1 - minBins + 1, x$region2 - minBins + 1)] <- x$p.adj
    mat[cbind(x$region2 - minBins + 1, x$region1 - minBins + 1)] <- x$p.adj
    allBreaks <- seq(-0.694, 0.3, by = 0.095)
    labelsPval <- round((exp(allBreaks) - 0.5), digits = 3)
    suppressMessages({
      pp <- plotSim(mat) + scale_fill_viridis(
        name = "adj p-value", breaks = allBreaks,
        labels = labelsPval, direction = -1
      )
    })
  }

  if (whichPlot == "logFC") {
    mat <- matrix(0, p, p)
    mat[cbind(x$region1 - minBins + 1, x$region2 - minBins + 1)] <- x$logFC
    mat[cbind(x$region2 - minBins + 1, x$region1 - minBins + 1)] <- x$logFC
    bInf <- min(mat)
    bSup <- max(mat)
    allBreaks <- seq(bInf, bSup, by = ((bSup - bInf)) / 8)
    labelsPval <- round(allBreaks, digits = 3)
    suppressMessages({
      pp <- plotSim(mat, log = FALSE) + scale_fill_viridis(
        name = "logFC", breaks = allBreaks,
        labels = labelsPval, direction = -1
      )
    })
  }
  return(pp)
}
