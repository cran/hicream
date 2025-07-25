## ----loadLib------------------------------------------------------------------
library("hicream")
# checking if Python and Python modules are available to avoid errors in vignette building
modules_avail <- reticulate::py_available(initialize = TRUE) &&
  reticulate::py_module_available("sklearn") &&
  reticulate::py_module_available("kneebow") &&
  reticulate::py_module_available("pandas") &&
  reticulate::py_module_available("numpy")

## ----loadData-----------------------------------------------------------------
replicates <- 1:2
cond <- "90"
allBegins <- interaction(expand.grid(replicates, cond), sep = "-")
allBegins <- as.character(allBegins)
chromosome <- 1
nbChr <- 1
allMat <- sapply(allBegins, function(ab) {
  matFile <- paste0("Rep", ab, "-chr", chromosome, "_200000.bed")
})
index <- system.file("extdata", "index.200000.longest18chr.abs.bed",
                     package = "hicream")
format <- rep("HiC-Pro", length(replicates) * length(cond) * nbChr)
binsize <- 200000
files <- system.file("extdata", unlist(allMat), package = "hicream")
exData <- loadData(files, index, chromosome, normalize = TRUE)

## ----pighic-------------------------------------------------------------------
data("pighic")
head(pighic)

## ----performTest, hold=TRUE---------------------------------------------------
resdiff <- performTest(pighic$data, pighic$conditions)
resdiff
summary(resdiff)

## ----plotPerformTest----------------------------------------------------------
plot(resdiff)
plot(resdiff, which_plot = "p.adj")
plot(resdiff, which_plot = "logFC")

## ----AggloClust2D-------------------------------------------------------------
if (modules_avail) {
  res2D <- AggloClust2D(pighic$data)
  res2D 
  summary(res2D)
}

## ----plotAggloClust2D---------------------------------------------------------
if (modules_avail) {
  plot(res2D)
}

## ----postHoc------------------------------------------------------------------
if (modules_avail) {
  clusters <- res2D$clustering
  alpha <- 0.05
  resposthoc <- postHoc(resdiff, clusters, alpha)
  resposthoc
  summary(resposthoc)
}

## ----plotPostHoc--------------------------------------------------------------
if (modules_avail) plot(resposthoc)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

