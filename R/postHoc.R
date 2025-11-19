#' @title Perform Post Hoc Inference
#'
#' @description This function performs post-hoc inference on all provided clusters
#' using \code{\link{performTest}} results.
#'
#' @param resdiff An object of class \code{resdiff} obtained from the
#' function \code{\link{performTest}}.
#' @param clusters A vector corresponding to a clustering of resdiff rows.
#' @param alpha A number between 0 and 1 at which the computed post hoc bounds
#' will be valid.
#' @param fill A boolean value. If TRUE (the default), enforce the total 
#' number of tests to be \emph{p*(p+1)/2} by adding ones to the p-value vector 
#' for non tested entries.
#' @param stepDown A boolean value. If FALSE, the single step Simes method is
#' used. If TRUE (the default), the step-down Simes method is used. The latter
#' is more powerful when the signal is strong.
#'
#' @return An object of class \code{resposthoc} containing a matrix with true
#' positive proportions for each interaction and a dataframe with the following
#' entries:
#' \item{region1}{The first bin of the interaction.}
#' \item{region2}{The second bin of the interaction.}
#' \item{clust}{The cluster the interaction belongs to.}
#' \item{TPRate}{The minimal post hoc true positive proportion of the cluster
#' the interaction belongs to.}
#' \item{p.value}{The p-value of the diffHic test.}
#' \item{p.adj}{The adjusted p-value of the diffHic test.}
#' \item{logFC}{The log2-fold-change of the interaction.}
#' \item{meanlogFC}{The mean of the log2-fold-change for the cluster the
#' interaction belongs to.}
#' \item{varlogFC}{The variance of the log2-fold-change for the cluster the
#' interaction belongs to.}
#' \item{propPoslogFC}{The proportion of interactions with positive
#' log2-fold-change in the cluster the interaction belongs to.}
#'
#' @author Élise Jorge \email{elise.jorge@inrae.fr}\cr
#' Sylvain Foissac \email{sylvain.foissac@inrae.fr}\cr
#' Toby Dylan Hocking \email{toby.hocking@r-project.org}\cr
#' Pierre Neuvial \email{pierre.neuvial@math.univ-toulouse.fr}\cr
#' Nathalie Vialaneix \email{nathalie.vialaneix@inrae.fr}
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom dplyr ungroup
#' @importFrom limma squeezeVar
#' @importFrom reshape2 colsplit
#' @importFrom stats cophenetic
#' @importFrom stats pt
#' @importFrom rlang .data
#' @importFrom viridis scale_fill_viridis
#' @importFrom utils write.csv
#' @importFrom stats var
#' @importFrom utils head
#' @importFrom rlang .data
#'
#' @examples
#' data("pighic")
#' resdiff <- performTest(pighic$data, pighic$conditions)
#' \donttest{res2D <- AggloClust2D(pighic$data)
#' if (!is.null(res2D)) { # in case Python or modules are not available
#'   clusters <- res2D$clustering
#'   alpha <- 0.05
#'   resposthoc <- postHoc(resdiff, clusters, alpha)
#'   resposthoc
#'   summary(resposthoc)
#'   plot(resposthoc)
#' }}

postHoc <- function(resdiff, clusters, alpha, fill = FALSE, stepDown = TRUE) {
  uniqueBins <- unique(c(resdiff$region1, resdiff$region2))
  p <- max(uniqueBins) - min(uniqueBins) + 1
  matPostHoc <- matrix(0, p, p)
  diffPvalInd <- c(resdiff$p.value, 1)
  matPval <- matrix(length(diffPvalInd), p, p) # Untested interactions
  # have a p-value of one. Put at last position in p-value vector.
  # matPval[i,j] is the position of (i,j) interaction in diffHic vector.
  matPval[cbind(
    resdiff$region1 - min(uniqueBins) + 1,
    resdiff$region2 - min(uniqueBins) + 1
  )] <- seq(1:length(resdiff$region1))
  matPval[cbind(
    resdiff$region2 - min(uniqueBins) + 1,
    resdiff$region1 - min(uniqueBins) + 1
  )] <- seq(1:length(resdiff$region1))
  vecRatioVPClust <- c()
  ## Put clustering in list format
  pixels <- seq(1, length(clusters))
  clusters <- sapply(seq(1:length(unique(clusters))), function(ab) {
    elem <- pixels[clusters == ab]
    return(elem)
  })
  ## define Simes thresholds
  pvals <- resdiff$p.value
  if (fill) {  # add 1:s for untested interactions
    m <- p * (p + 1) / 2 # Total number of interactions to consider in theory
    nb_ones <- m - length(resdiff$p.value)
    pvals <- c(pvals, rep(1, nb_ones)) 
  }
  thrSimes <- getSimesThresholds(pvals, alpha, stepDown = stepDown)
  for (i in 1:length(clusters)) {
    clusterElem <- clusters[[i]]
    clustBins1 <- resdiff$region1[clusterElem] - min(uniqueBins) + 1
    clustBins2 <- resdiff$region2[clusterElem] - min(uniqueBins) + 1
    indPval <- matPval[cbind(clustBins1, clustBins2)]
    # Compute upper bound on number of false positives
    FP <- curveMaxFP(diffPvalInd[indPval], thrSimes)
    tot <- length(clusterElem)
    # Ratio of True Positives
    ratio <- 1 - (FP[length(FP)] / tot)
    matPostHoc[cbind(clustBins1, clustBins2)] <- ratio
    matPostHoc[cbind(clustBins2, clustBins1)] <- ratio
    vecRatioVPClust[i] <- ratio
  }

  ## Create output table with cluster information
  lenClust <- lengths(clusters)
  clust <- unlist(lapply(c(1:length(clusters)), function(i) {
    return(rep(i, lenClust[i]))
  }))
  postHocRes <- unlist(lapply(c(1:length(clusters)), function(i) {
    return(rep(vecRatioVPClust[i], lenClust[i]))
  }))
  dfRes <- data.frame("int" = unlist(clusters))
  dfRes$clust <- clust
  dfRes$TPRate <- postHocRes
  resdiff <- do.call(cbind.data.frame, resdiff)
  resdiff$X <- seq(1:length(resdiff$region1))
  dfRes <- merge(dfRes, resdiff, by.x = "int", by.y = "X")
  dfRes <- dfRes[, c(4, 5, 2, 3, 6, 7, 8)]

  ## LogFC
  ## Compute logFC mean, variance and sign.
  ## Add mean column.
  dfRes <- dfRes %>%
    group_by(clust) %>%
    mutate(meanlogFC = mean(.data$logFC))

  ## Add variance column.
  dfRes <- dfRes %>%
    group_by(clust) %>%
    mutate(varlogFC = var(.data$logFC))

  ## Add proportion of positive logFC column.
  dfRes <- dfRes %>%
    group_by(clust) %>%
    mutate(propPoslogFC = ((sum(.data$logFC >= 0)) / (length(.data$logFC))))

  dfRes <- dfRes[with(dfRes, order(region1, region2)), ]

  outList <- list(matPostHoc, dfRes)
  names(outList) <- c("mat", "dfres")
  # Assign class to the list
  class(outList) <- "resposthoc"

  return(outList)
}



##### Methods for resposthoc object ####
#' @exportS3Method
#' @param x a \code{resposthoc} object to print
#' @param ... not used
#' @rdname postHoc


print.resposthoc <- function(x, ...) {
  # print title
  cat("Posthoc results.")
  cat("\n\n")

  ## Print the dataframe
  print(x$dfres, n = 10)
}


#' @exportS3Method
#' @param object a \code{resposthoc} object to summarize
#' @param ... not used
#' @rdname postHoc

summary.resposthoc <- function(object, ...) {
  # print title
  cat("Summary of post hoc results.")
  cat("\n\n")

  ## Only look at TPRate, p.vlaue, p.adj, logFC and propPoslogFC.
  ## Print the dataframe
  dfSum <- object$dfres %>% ungroup()
  dfSum <- dfSum %>%
    select(
      .data$TPRate, .data$p.value, .data$p.adj, .data$logFC,
      .data$propPoslogFC
    )
  summary(dfSum)
}


#' @exportS3Method
#' @param x a \code{resposthoc} object to plot
#' @param ... not used
#' @rdname postHoc


plot.resposthoc <- function(x, ...) {
  matposthoc <- x$mat
  maxRatio <- max(matposthoc)
  allBreaks <- seq(0, maxRatio, by = round(maxRatio / 6, digits = 2))
  labels <- round((exp(allBreaks) - 0.5) * 100)

  # Plot post-hoc ratios
  suppressMessages({
  phRatios <- plotSim(matposthoc, log = F) +
    scale_fill_viridis(
      name = "TDP",
      breaks = allBreaks,
      labels = allBreaks,
      direction = 1
    )
  })
  return(phRatios)
}
