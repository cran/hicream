##### for loadData
#' @importFrom GenomeInfoDb seqlevels
mergeInteractionSet <- function(interactionSet1, interactionSet2, fill = 0) {
  unionInteractions <- GenomicRanges::union(
    interactions(interactionSet1),
    interactions(interactionSet2)
  )
  interactionSet1 <- fillInteractionSet(
    interactionSet1,
    unionInteractions,
    fill
  )
  interactionSet2 <- fillInteractionSet(
    interactionSet2,
    unionInteractions,
    fill
  )
  newiset <- BiocGenerics::cbind(interactionSet1, interactionSet2)
  return(newiset)
}

#' @importFrom utils str
fillInteractionSet <- function(interactionSet, interactionSetUnion, fill = 0) {
  a <- seqlevels(interactionSetUnion)
  b <- seqlevels(interactionSet)
  message(str(union(setdiff(a, b), setdiff(b, a))))
  seqlevels(interactionSet) <- unique(c(
    seqlevels(interactionSetUnion),
    seqlevels(interactionSet)
  ))
  seqlevels(interactionSetUnion) <- seqlevels(interactionSet)
  over <- GenomicRanges::match(interactionSet, interactionSetUnion)
  totalColumns <- ncol(interactionSet)
  newAssays <- matrix(
    rep(fill, length(interactionSetUnion) * totalColumns),
    ncol = totalColumns
  )
  newAssays[over, ] <- SummarizedExperiment::assay(interactionSet)
  out <- InteractionSet::InteractionSet(
    newAssays,
    interactionSetUnion,
    colData = SummarizedExperiment::colData(interactionSet)
  )
  return(out)
}


#### for AggloClust2d ####
fromDGE2Counts <- function(curDge, offset = NULL) {
  counts <- curDge@assays@data@listData$counts
  libSizes <- colSums(counts)
  if (is.null(offset)) {
    if ("offset" %in% names(curDge@assays@data@listData)) {
      offsets <- curDge@assays@data@listData$offset
      offsets <- offsets - mean(log(libSizes))
    } else {
      offsets <- -mean(log(libSizes))
    }
  }
  counts <- counts / exp(offsets)

  countMatrix <- data.frame(interactions(curDge))
  countMatrix <- data.frame(
    "bin1" = countMatrix$anchor1.id,
    "bin2" = countMatrix$anchor2.id,
    counts
  )

  return(countMatrix)
}

neighborsToMat <- function(neighbors, nbInt) {
  return(sparseMatrix(
    i = c(neighbors[, 1], neighbors[, 2]),
    j = c(neighbors[, 2], neighbors[, 1]),
    x = neighbors[, 3], dims = c(nbInt, nbInt)
  ))
}

#' @importFrom auk auk_get_awk_path
#' @importFrom utils write.table read.table
buildNeighborsR <- function(counts) {
  newDf <- data.frame(bin1 = counts$bin1, bin2 = counts$bin2)
  inFile <- tempfile(fileext = ".bed")
  write.table(newDf,
    file = inFile, row.names = FALSE, col.names = FALSE,
    sep = "\t", quote = FALSE
  )
  outFile <- tempfile(fileext = ".bed")

  awkPath <- auk_get_awk_path()
  if (is.na(awkPath)) {
    stop("'buildNeighborsR2' requires a valid AWK install.")
  }

  # run command
  awkClean <- "'left[$1]{print left[$1],NR,1}up[$2]{print up[$2],NR,1}{left[$1]=NR;up[$2]=NR}' OFS=\"\t\""
  exitCode <- system2(awkPath, paste(awkClean, inFile), outFile)
  if (exitCode == 0) {
    outFile
  } else {
    exitCode
    stop("stopping here... error encountered in AWK command")
  }

  neighbors <- read.table(outFile)
  file.remove(inFile)
  file.remove(outFile)
  return(neighbors)
}

fromAgglomerative2hclust <- function(merge, height, call) {
  height <- height[, 1]
  hierarchy <- merge + 1 # +1 fixes index shift b/w R and python
  hierMat <- as.matrix(hierarchy)
  dimnames(hierMat) <- NULL
  nbInt <- nrow(hierMat) + 1
  leafInd <- hierMat <= nbInt
  clustInd <- hierMat > nbInt
  hierMat[leafInd] <- -hierMat[leafInd]
  hierMat[clustInd] <- hierMat[clustInd] - nbInt
  hc <- list(
    merge = hierMat,
    height = height,
    order = seq(length = nbInt),
    labels = seq(length = nbInt),
    method = "Constrained HC with Ward linkage from sklearn"
  )
  mergeVect <- as.vector(t(hc[["merge"]]))
  hc[["order"]] <- -mergeVect[mergeVect < 0]
  hc$call <- call
  hc$"dist.method" <- "euclidean"
  class(hc) <- "hclust"
  return(hc)
}

#### for postHoc ####

# Upper bound for the number of false discoveries among most significant items
# Taken from R package 'sousoucis' https://github.com/sanssouci-org/sanssouci
curveMaxFP <- function(p.values, thr) {
  s <- length(p.values)
  if (s == 0) {
    return(numeric(0L))
  }
  p.values <- sort(p.values)
  thr <- sort(thr)
  
  kMax <- length(thr)
  if (s < kMax){  # truncate thr to first 's' values
    seqK <- seq(from = 1, to = s, by = 1)
    thr <- thr[seqK]
  } else { # complete 'thr' to length 's' with its last value
    thr <- c(thr, rep(thr[kMax], s - kMax))
  }
  ## sanity checks
  stopifnot(length(thr) == s)
  rm(kMax)
  
  K <- rep(s, s) ## K[i] = number of k/ T[i] <= s[k]
  Z <- rep(s, s) ## Z[k] = number of i/ T[i] >  s[k] = cardinal of R_k
  ## 'K' and 'Z' are initialized to their largest possible value (both 's')
  kk <- 1
  ii <- 1
  while ((kk <= s) && (ii <= s)) {
    if (thr[kk] > p.values[ii]) {
      K[ii] <- kk-1
      ii <- ii+1
    } else {
      Z[kk] <- ii-1
      kk <- kk+1
    }
  }
  Vbar <- numeric(s)
  ww <- which(K > 0)
  A <- Z - (1:s) + 1
  cA <- cummax(A)[K[ww]]  # cA[i] = max_{k<K[i]} A[k]
  Vbar[ww] <- pmin(ww - cA, K[ww])
  
  return(Vbar)
}

# Upper bound for the number of false null hypotheses among most significant 
# items
# Note: roxygen comments not interpreted but kept for dev information
#
# @param p.values A vector containing m p-values
# @param alpha A numeric value between 0 and 1, the targed risk level
# @param stepDown A boolean value. If TRUE (the default), the adaptive
# (step-down) method is used. If FALSE, the single-step method is used.
# @return thr A vector of \eqn{m_0} JER-controlling thresholds based on the Simes
# family of the form \code{alpha*k/m_0} for \eqn{1 \leq k \leq m_0}
# @details The default adaptive (step-down) method is at least as powerful (and
# possibly more) than the single-step method. Its complexity is still $O(m)$
# since it applies 'curveMaxFP' a small number of times. Therefore it should
# be preferred.
# @examples
#stat <- c(rnorm(1000, mean = c(1:1000)/1000*5, sd = 1), 
#          rnorm(100, mean = 0, sd = 1))
# pvals <- 2*(1 - pnorm(abs(stat)))
# thr <- hicream:::getSimesThresholds(pvals, 0.05)
# str(thr)
# tail(hicream:::curveMaxFP(pvals, thr))
# thr <- hicream:::getSimesThresholds(pvals, 0.05, stepDown = FALSE)
# str(thr)
# tail(hicream:::curveMaxFP(pvals, thr))}

getSimesThresholds <- function(p.values, alpha, stepDown = TRUE) {
  m <- length(p.values)
  thrSimes <- alpha * (1:m) / m # single-step
  if (stepDown) {
    m0_hat0 <- m
    FP <- curveMaxFP(p.values, thrSimes)
    m0_hat <- FP[m]
    thrSimes <- alpha * (1:m0_hat) / m0_hat # one step down
    while (m0_hat < m0_hat0) { # stepping down as needed
      m0_hat0 <- m0_hat
      FP <- curveMaxFP(p.values, thrSimes)
      m0_hat <- FP[m]
      thrSimes <- alpha * (1:m0_hat) / m0_hat
    }
  }
  
  return(thrSimes)
}

