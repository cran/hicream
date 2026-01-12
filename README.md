---
title: "Introduction to hicream"
author: "Elise Jorge, Toby Hocking, Sylvain Foissac, Pierre Neuvial, and Nathalie Vialaneix"
date: "2026-01-10"
output: 
  html_document:
    toc: yes
    keep_md: TRUE
    self_contained: yes
vignette: >
  %\VignetteIndexEntry{Introduction to hicream}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!-- badges: start -->
[![CRAN version](https://www.r-pkg.org/badges/version/hicream)](https://cran.r-project.org/package=hicream)
[![CRAN checks](https://badges.cranchecks.info/summary/hicream.svg)](https://cran.r-project.org/web/checks/check_results_hicream.html)
[![CRANLOGS](https://cranlogs.r-pkg.org/badges/hicream)](https://cran.r-project.org/package=hicream)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:c129fb55ad150cf47c1b2f247d5e3de7d250b58d/)](https://archive.softwareheritage.org/swh:1:dir:c129fb55ad150cf47c1b2f247d5e3de7d250b58d;origin=https://forge.inrae.fr/scales/hicream)
<!-- badges: end -->

# Introduction 

`hicream` is an **R** package designed to perform Hi-C data differential 
analysis. More specifically, it performs a pixel-level differential analysis 
using `diffHic` and, using a two-dimensional connectivity constrained 
clustering, renders a post hoc bound on True Discovery Proportion for each 
cluster. This method allows to identify differential genomic regions and 
quantifies signal in those regions. 


``` r
library("hicream")
```

```
## Loading required package: reticulate
```

```
## Warning: multiple methods tables found for 'seqinfo'
```

```
## Warning: multiple methods tables found for 'seqinfo<-'
```

```
## Warning: multiple methods tables found for 'seqnames'
```

```
## Warning: multiple methods tables found for 'seqnames<-'
```

```
## Warning: multiple methods tables found for 'commonName'
```

```
## Warning: multiple methods tables found for 'provider'
```

```
## Warning: multiple methods tables found for 'providerVersion'
```

```
## Warning: multiple methods tables found for 'releaseDate'
```


## How to load external data?

Using `hicream` requires to load data that corresponds to Hi-C matrices and 
their index, the latter in bed format. In the following code, the paths to Hi-C 
matrices and the index file path are used in the `loadData` function alongside
the chromosome number. The option `normalize = TRUE` allows to perform a cyclic 
LOESS normalization. The output is an `InteractionSet` object. 


``` r
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
```

```
##  chr(0)
```

```
## 
```

```
##  chr(0)
```

```
## 
```


## `pighic` dataset 

The `pighic` dataset has been produced using 6 Hi-C matrices (3 in each 
condition) obtained from two different developmental stages of pig embryos (90 
and 110 days of gestation) corresponding to chromosome 1. The dataset contains 
two elements, the first is an output of the `loadData` function corresponding to
normalized data. The second is the vector giving, for each matrix, the condition
it belongs to. 


``` r
data("pighic")
head(pighic)
```

```
## $data
## class: InteractionSet 
## dim: 21 6 
## metadata(1): width
## assays(2): counts offset
## rownames: NULL
## rowData names(2): anchor1.id anchor2.id
## colnames: NULL
## colData names(1): totals
## type: ReverseStrictGInteractions
## regions: 6
## 
## $conditions
## [1] "90"  "90"  "90"  "110" "110" "110"
```


# Perform `hicream` test

Once data have been loaded, the three steps of the analysis can be performed. 

## Use `performTest` function

In this first step, the output of a `loadData` object is used, with a vector 
indicating the condition of each sample. Here, we use data stored in `pighic`. 
`performTest` outputs the result of the `diffHic` pixel-level differential 
analysis. 


``` r
resdiff <- performTest(pighic$data, pighic$conditions)
resdiff
```

```
## Testing for differential interactions using diffHic.
## 
##   region1 region2   p.value     p.adj       logFC
## 1    1125    1125 0.7318482 0.8538229 -0.01830721
## 2    1125    1126 0.0562111 0.2342558 -0.16338780
##  [ reached 'max' / getOption("max.print") -- omitted 19 rows ]
```

``` r
summary(resdiff)
```

```
## Summary of diffHic results.
```

```
##     p.value              p.adj              logFC          
##  Min.   :0.0002064   Min.   :0.004334   Min.   :-0.163388  
##  1st Qu.:0.0669302   1st Qu.:0.234256   1st Qu.:-0.034500  
##  Median :0.2751093   Median :0.525209   Median : 0.003099  
##  Mean   :0.3960997   Mean   :0.580792   Mean   : 0.022271  
##  3rd Qu.:0.6924986   3rd Qu.:0.853823   3rd Qu.: 0.040991  
##  Max.   :0.9537715   Max.   :0.953771   Max.   : 0.349596
```

Several plotting options allow to look at the $p$-value map, adjusted $p$-value
map or log-fold-change map. 


``` r
plot(resdiff)
```

![](README_files/figure-html/plotPerformTest-1.png)<!-- -->

``` r
plot(resdiff, whichPlot = "p.adj")
```

![](README_files/figure-html/plotPerformTest-2.png)<!-- -->

``` r
plot(resdiff, whichPlot = "logFC")
```

![](README_files/figure-html/plotPerformTest-3.png)<!-- -->


## Perform two-dimensional connectivity-constrained clustering

The `AggloClust2D` function uses either the output of `loadData` 
or `resdiff` function. If the output of `loadData` is used, the clustering is 
performed on the counts. Using the output of `resdiff` amounts to 
performing clustering on adjusted $p$-values and log-fold-changes. 
In both cases, a connectivity graph of pixels is built and a two-dimensional
hierarchical clustering under the connectivity constraint is performed. 
The function renders a clustering 
corresponding to the optimal number of clusters found by the elbow heuristic. 
However, a clustering corresponding to any other number of clusters (chosen by
the user) can be obtained by specifying a value for the input argument 
`nbClust`. 


``` r
res2D_counts <- AggloClust2D(pighic$data)
res2D_counts 
```

```
## Tree obtained from constrained 2D clustering.
## 
## Call:
## AggloClust2D(pighic$data)
## 
## Cluster method   : Constrained HC with Ward linkage from sklearn 
## Distance         : euclidean 
## Number of objects: 21 
## 
## 
## Optimal number of clusters: 7 
## 
## Clustering:
## 1 1 5 2 2 2 1 1 2 4 ...
```

``` r
summary(res2D_counts)
```

```
## Summary of 2D constrained clustering results.
```

```
##             Length Class  Mode     
## merge       40     -none- numeric  
## height      20     -none- numeric  
## order       21     -none- numeric  
## labels      21     -none- numeric  
## method       1     -none- character
## call         2     -none- call     
## dist.method  1     -none- character
```


``` r
res2D_diff <- AggloClust2D(resdiff)
res2D_diff 
```

```
## Tree obtained from constrained 2D clustering.
## 
## Call:
## AggloClust2D(resdiff)
## 
## Cluster method   : Constrained HC with Ward linkage from sklearn 
## Distance         : euclidean 
## Number of objects: 21 
## 
## 
## Optimal number of clusters: 3 
## 
## Clustering:
## 1 1 1 1 3 1 1 1 1 1 ...
```

``` r
summary(res2D_diff)
```

```
## Summary of 2D constrained clustering results.
```

```
##             Length Class  Mode     
## merge       40     -none- numeric  
## height      20     -none- numeric  
## order       21     -none- numeric  
## labels      21     -none- numeric  
## method       1     -none- character
## call         2     -none- call     
## dist.method  1     -none- character
```

Plotting the output shows the tree that corresponds to the hierarchy of 
clusters.


``` r
plot(res2D_diff)
```

![](README_files/figure-html/plotAggloClust2D-1.png)<!-- -->


## Perform post hoc inference

Post hoc inference is performed by the `postHoc` function, using the output 
of `performTest`, and a clustering (the latter can be obtained with 
`AggloClust2D`). The user sets a level of test confidence $\alpha$, typically 
equal to $0.05$. 


``` r
clusters <- res2D_diff$clustering
alpha <- 0.05
resposthoc <- postHoc(resdiff, clusters, alpha)
resposthoc
```

```
## Posthoc results.
## 
## # A tibble: 21 × 10
## # Groups:   clust [3]
##    region1 region2 clust TPRate p.value p.adj    logFC meanlogFC varlogFC
##      <int>   <int> <int>  <dbl>   <dbl> <dbl>    <dbl>     <dbl>    <dbl>
##  1    1125    1125     1      0  0.732  0.854 -0.0183   -0.00828  0.00467
##  2    1125    1126     1      0  0.0562 0.234 -0.163    -0.00828  0.00467
##  3    1125    1127     1      0  0.800  0.884 -0.00890  -0.00828  0.00467
##  4    1125    1128     1      0  0.882  0.926 -0.00775  -0.00828  0.00467
##  5    1125    1129     3      0  0.0313 0.234  0.275     0.275   NA      
##  6    1125    1130     1      0  0.269  0.525 -0.0659   -0.00828  0.00467
##  7    1126    1126     1      0  0.692  0.854  0.0280   -0.00828  0.00467
##  8    1126    1127     1      0  0.186  0.489 -0.0648   -0.00828  0.00467
##  9    1126    1128     1      0  0.409  0.716  0.0364   -0.00828  0.00467
## 10    1126    1129     1      0  0.251  0.525  0.0487   -0.00828  0.00467
## # ℹ 11 more rows
## # ℹ 1 more variable: propPoslogFC <dbl>
```

``` r
summary(resposthoc)
```

```
## Summary of post hoc results.
```

```
##      TPRate           p.value              p.adj              logFC          
##  Min.   :0.00000   Min.   :0.0002064   Min.   :0.004334   Min.   :-0.163388  
##  1st Qu.:0.00000   1st Qu.:0.0669302   1st Qu.:0.234256   1st Qu.:-0.034500  
##  Median :0.00000   Median :0.2751093   Median :0.525209   Median : 0.003099  
##  Mean   :0.04762   Mean   :0.3960997   Mean   :0.580792   Mean   : 0.022271  
##  3rd Qu.:0.00000   3rd Qu.:0.6924986   3rd Qu.:0.853823   3rd Qu.: 0.040991  
##  Max.   :1.00000   Max.   :0.9537715   Max.   :0.953771   Max.   : 0.349596  
##   propPoslogFC   
##  Min.   :0.4737  
##  1st Qu.:0.4737  
##  Median :0.4737  
##  Mean   :0.5238  
##  3rd Qu.:0.4737  
##  Max.   :1.0000
```

Plotting the output of `postHoc` function shows the minimal proportion of True 
Discoveries in each cluster. 


``` r
plot(resposthoc)
```

![](README_files/figure-html/plotPostHoc-1.png)<!-- -->


# Session Information


``` r
sessionInfo()
```

```
## R version 4.5.2 (2025-10-31)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Paris
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] hicream_0.0.2     reticulate_1.44.1
## 
## loaded via a namespace (and not attached):
##   [1] bitops_1.0-9                gridExtra_2.3              
##   [3] rlang_1.1.6                 magrittr_2.0.4             
##   [5] otel_0.2.0                  matrixStats_1.5.0          
##   [7] compiler_4.5.2              png_0.1-8                  
##   [9] vctrs_0.6.5                 reshape2_1.4.5             
##  [11] stringr_1.6.0               pkgconfig_2.0.3            
##  [13] crayon_1.5.3                fastmap_1.2.0              
##  [15] XVector_0.50.0              labeling_0.4.3             
##  [17] Rsamtools_2.26.0            rmarkdown_2.30             
##  [19] UCSC.utils_1.6.1            xfun_0.55                  
##  [21] cachem_1.1.0                cigarillo_1.0.0            
##  [23] GenomeInfoDb_1.46.2         jsonlite_2.0.0             
##  [25] rhdf5filters_1.22.0         DelayedArray_0.36.0        
##  [27] Rhdf5lib_1.32.0             BiocParallel_1.44.0        
##  [29] parallel_4.5.2              R6_2.6.1                   
##  [31] bslib_0.9.0                 stringi_1.8.7              
##  [33] RColorBrewer_1.1-3          limma_3.66.0               
##  [35] rtracklayer_1.70.0          GenomicRanges_1.62.1       
##  [37] jquerylib_0.1.4             Rcpp_1.1.0                 
##  [39] Seqinfo_1.0.0               SummarizedExperiment_1.40.0
##  [41] knitr_1.51                  IRanges_2.44.0             
##  [43] Matrix_1.7-4                splines_4.5.2              
##  [45] tidyselect_1.2.1            rstudioapi_0.17.1          
##  [47] abind_1.4-8                 yaml_2.3.12                
##  [49] viridis_0.6.5               codetools_0.2-20           
##  [51] curl_7.0.0                  lattice_0.22-7             
##  [53] tibble_3.3.0                plyr_1.8.9                 
##  [55] InteractionSet_1.38.0       Biobase_2.70.0             
##  [57] withr_3.0.2                 S7_0.2.1                   
##  [59] csaw_1.44.0                 evaluate_1.0.5             
##  [61] capushe_1.1.3               Biostrings_2.78.0          
##  [63] pillar_1.11.1               MatrixGenerics_1.22.0      
##  [65] auk_0.9.0                   stats4_4.5.2               
##  [67] generics_0.1.4              rprojroot_2.1.1            
##  [69] RCurl_1.98-1.17             S4Vectors_0.48.0           
##  [71] ggplot2_4.0.1               sparseMatrixStats_1.22.0   
##  [73] scales_1.4.0                glue_1.8.0                 
##  [75] metapod_1.18.0              tools_4.5.2                
##  [77] dendextend_1.19.1           BiocIO_1.20.0              
##  [79] BSgenome_1.76.0             locfit_1.5-9.12            
##  [81] GenomicAlignments_1.46.0    XML_3.99-0.20              
##  [83] rhdf5_2.54.1                grid_4.5.2                 
##  [85] edgeR_4.8.1                 adjclust_0.6.11            
##  [87] Rhtslib_3.6.0               restfulr_0.0.16            
##  [89] cli_3.6.5                   rappdirs_0.3.3             
##  [91] S4Arrays_1.10.1             viridisLite_0.4.2          
##  [93] dplyr_1.1.4                 gtable_0.3.6               
##  [95] sass_0.4.10                 digest_0.6.39              
##  [97] BiocGenerics_0.56.0         SparseArray_1.10.8         
##  [99] diffHic_1.40.0              rjson_0.2.23               
## [101] farver_2.1.2                htmltools_0.5.9            
## [103] lifecycle_1.0.4             httr_1.4.7                 
## [105] here_1.0.2                  statmod_1.5.1              
## [107] MASS_7.3-65
```
