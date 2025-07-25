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

# Upper bound for the number of false discoveries among most significant items
# Taken from R package 'sousoucis' (on github)
# Courtesy of Gilles Blanchard, Nicolas Enjalbert-Courrech, Pierre Neuvial, and
#   Etienne Roquain
# Reference: Enjalbert-Courrech, N. & Neuvial, P. (2022). Powerful and 
#   interpretable control of false discoveries in two-group differential 
#   expression studies. Bioinformatics. doi: 10.1093/bioinformatcs/btac693
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
  Vbar
}