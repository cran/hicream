#' @title Load and Normalize Hi-C Data
#'
#' @description This function loads data necessary for the analysis and outputs
#' them in a suitable format for \code{performTest} and \code{2Dclust}.
#'
#' @param files character vector. Paths to Hi-C matrices in bed format.
#' @param index character. A path to an index file in bed format.
#' @param chromosome character or integer. Chromosome to select.
#' @param normalize logical. Whether or not to normalize the output (with MA
#' method). Set to \code{TRUE} by default.
#'
#' @return An \code{\link[InteractionSet]{InteractionSet}} corresponding to all
#' interactions present in at least one of the input matrices and corresponding
#' counts across all matrices.
#'
#' @author Élise Jorge \email{elise.jorge@inrae.fr}\cr
#' Sylvain Foissac \email{sylvain.foissac@inrae.fr}\cr
#' Pierre Neuvial \email{pierre.neuvial@math.univ-toulouse.fr}\cr
#' Nathalie Vialaneix \email{nathalie.vialaneix@inrae.fr}
#'
#' @export
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame match
#' @importFrom InteractionSet GInteractions interactions InteractionSet
#' @importFrom InteractionSet `interactions<-`
#' @importFrom GenomeInfoDb `seqlevels<-`
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges union
#' @importFrom BiocGenerics cbind
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom csaw normOffsets
#' @importFrom methods as
#'
#' @examples
#' replicates <- 1:2
#' cond <- "90"
#' allBegins <- interaction(expand.grid(replicates, cond), sep = "-")
#' allBegins <- as.character(allBegins)
#' chromosome <- 1
#' nbChr <- 1
#' allMat <- sapply(allBegins, function(ab) {
#'   matFile <- paste0("Rep", ab, "-chr", chromosome, "_200000.bed")
#'   })
#' index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
#'                     package = "hicream")
#'                     format <- rep("HiC-Pro", length(replicates) * length(cond) * nbChr)
#' binsize <- 200000
#' files <- system.file("extdata", unlist(allMat), package = "hicream")
#' exData <- loadData(files, index, chromosome, normalize = TRUE)

loadData <- function(files, index, chromosome, normalize = TRUE) {
  # call HiCDOCDataSet function
  index <- read.table(index)
  names(index) <- c("chr", "start", "end", "id")
  # filter chromosome
  index <- index[index$chr == chromosome, ]
  chr <- unique(index$chr)
  resolution <- index[1, 3] - index[1, 2] + 1

  matrices <- lapply(files, function(filename) {
    amatrix <- read.table(filename)
    names(amatrix) <- c("bin1", "bin2", "count")
    bins1 <- amatrix[, 1]
    sel1 <- match(bins1, index$id)
    bins2 <- amatrix[, 2]
    sel2 <- match(bins2, index$id)
    selected <- !is.na(sel1) & !is.na(sel2)

    bins1 <- index[sel1[selected], ]
    bins1 <- makeGRangesFromDataFrame(bins1, keep.extra.columns = TRUE)
    seqlevels(bins1) <- as.character(chr)
    bins2 <- index[sel2[selected], ]
    bins2 <- makeGRangesFromDataFrame(bins2, keep.extra.columns = TRUE)
    seqlevels(bins2) <- as.character(chr)

    gi <- GInteractions(anchor1 = bins1, anchor2 = bins2)
    libSizes <- DataFrame(totals = sum(amatrix[, 3]))

    omatrix <- amatrix[selected, 3]
    omatrix <- matrix(data = omatrix, ncol = 1)
    out <- InteractionSet(list(counts = omatrix), gi, colData = libSizes)
    return(out)
  })
  matrices <- Reduce(mergeInteractionSet, matrices)

  names(matrices@assays) <- "counts"
  metadata(matrices)$width <- resolution
  interactions(matrices) <- as(
    interactions(matrices),
    "ReverseStrictGInteractions"
  )

  if (normalize) matrices <- normOffsets(matrices, se.out = TRUE)
  return(matrices)
}
