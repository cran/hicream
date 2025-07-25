---
title: "Introduction to hicream"
author: "Elise Jorge, Sylvain Foissac, Pierre Neuvial, et Nathalie Vialaneix"
date: "2025-04-02"
output: 
  html_document:
    toc: yes
    keep_md: yes
vignette: >
  %\VignetteIndexEntry{Introduction to hicream}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


# Introduction 

`hicream` is an **R** package designed to perform Hi-c data differential 
analysis. More specifically, it performs a pixel-level differential analysis 
using `diffHic` and, using a two-dimensionnal connectivity constrained 
clustering, renders a post hoc bound on True Discovery Proportion for each 
cluster. This method allows to identify differential genomic regions and 
quantifies signal in those regions. 


``` r
library("hicream")
```

```
## Loading required package: reticulate
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
##  3rd Qu.:0.6924986   3rd Qu.:0.853823   3rd Qu.: 0.040990  
##  Max.   :0.9537715   Max.   :0.953771   Max.   : 0.349596
```

Several plotting options allow to look at the $p$-value map, adjusted $p$-value
map or log-fold-change map. 


``` r
plot(resdiff)
```

![](README_files/figure-html/plotPerformTest-1.png)<!-- -->

``` r
plot(resdiff, which_plot = "p.adj")
```

![](README_files/figure-html/plotPerformTest-2.png)<!-- -->

``` r
plot(resdiff, which_plot = "logFC")
```

![](README_files/figure-html/plotPerformTest-3.png)<!-- -->


## Perform two-dimensionnal connectivity-constrained clustering

The `AggloClust2D` function uses the output of `loadData` function in order to 
build a connectivity graph of pixels and perform a two-dimensionnal hierarchical
clustering under the connectivity constraint. The function renders a clustering 
corresponding to the optimal number of clusters found by the elbow heuristic. 
However, a clustering corresponding to any other number of clusters (chosen by
the user) can be obtained by specifying a value for the input argument 
`nbClust`. 


``` r
res2D <- AggloClust2D(pighic$data)
res2D 
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
summary(res2D)
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
plot(res2D)
```

![](README_files/figure-html/plotAggloClust2D-1.png)<!-- -->


## Perform post hoc inference

Post hoc inference is performed by the `postHoc` function, using the output 
of `performTest`, and a clustering (the latter can be obtained with 
`AggloClust2D`). The user sets a level of test confidence $\alpha$, typically 
equal to $0.05$. 


``` r
clusters <- res2D$clustering
alpha <- 0.05
resposthoc <- postHoc(resdiff, clusters, alpha)
resposthoc
```

```
## Posthoc results.
```

```
## Warning: `...` must be empty in `format.tbl()`
## Caused by error in `format_tbl()`:
## ! `...` must be empty.
## ✖ Problematic argument:
## • max = 10
```

```
## # A tibble: 21 × 10
## # Groups:   clust [7]
##    region1 region2 clust TPRate p.value p.adj    logFC meanlogFC varlogFC
##      <int>   <int> <int>  <dbl>   <dbl> <dbl>    <dbl>     <dbl>    <dbl>
##  1    1125    1125     1  0      0.732  0.854 -0.0183   -0.0546   0.00670
##  2    1125    1126     1  0      0.0562 0.234 -0.163    -0.0546   0.00670
##  3    1125    1127     5  0      0.800  0.884 -0.00890  -0.00890 NA      
##  4    1125    1128     2  0      0.882  0.926 -0.00775   0.0595   0.0225 
##  5    1125    1129     2  0      0.0313 0.234  0.275     0.0595   0.0225 
##  6    1125    1130     2  0      0.269  0.525 -0.0659    0.0595   0.0225 
##  7    1126    1126     1  0      0.692  0.854  0.0280   -0.0546   0.00670
##  8    1126    1127     1  0      0.186  0.489 -0.0648   -0.0546   0.00670
##  9    1126    1128     2  0      0.409  0.716  0.0364    0.0595   0.0225 
## 10    1126    1129     4  0.167  0.251  0.525  0.0487    0.0447   0.0266 
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
##  3rd Qu.:0.16667   3rd Qu.:0.6924986   3rd Qu.:0.853823   3rd Qu.: 0.040990  
##  Max.   :0.16667   Max.   :0.9537715   Max.   :0.953771   Max.   : 0.349596  
##   propPoslogFC   
##  Min.   :0.0000  
##  1st Qu.:0.5000  
##  Median :0.5000  
##  Mean   :0.5238  
##  3rd Qu.:0.6667  
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
## R version 4.4.3 (2025-02-28)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.2 LTS
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
## [1] hicream_0.0.1       reticulate_1.41.0.1
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_1.2.1            viridisLite_0.4.2          
##   [3] farver_2.1.2                dplyr_1.1.4                
##   [5] viridis_0.6.5               Biostrings_2.74.1          
##   [7] bitops_1.0-9                RCurl_1.98-1.16            
##   [9] fastmap_1.2.0               GenomicAlignments_1.42.0   
##  [11] XML_3.99-0.18               digest_0.6.37              
##  [13] lifecycle_1.0.4             capushe_1.1.2              
##  [15] statmod_1.5.0               magrittr_2.0.3             
##  [17] compiler_4.4.3              rlang_1.1.5                
##  [19] sass_0.4.9                  tools_4.4.3                
##  [21] yaml_2.3.10                 rtracklayer_1.66.0         
##  [23] knitr_1.50                  Rhtslib_3.2.0              
##  [25] labeling_0.4.3              S4Arrays_1.6.0             
##  [27] curl_6.2.1                  DelayedArray_0.32.0        
##  [29] plyr_1.8.9                  abind_1.4-8                
##  [31] BiocParallel_1.40.0         withr_3.0.2                
##  [33] BiocGenerics_0.52.0         grid_4.4.3                 
##  [35] stats4_4.4.3                colorspace_2.1-1           
##  [37] Rhdf5lib_1.28.0             edgeR_4.4.2                
##  [39] ggplot2_3.5.1               scales_1.3.0               
##  [41] MASS_7.3-65                 SummarizedExperiment_1.36.0
##  [43] cli_3.6.4                   rmarkdown_2.29             
##  [45] crayon_1.5.3                auk_0.8.0                  
##  [47] generics_0.1.3              metapod_1.14.0             
##  [49] rstudioapi_0.17.1           reshape2_1.4.4             
##  [51] rjson_0.2.23                httr_1.4.7                 
##  [53] rhdf5_2.50.2                cachem_1.1.0               
##  [55] stringr_1.5.1               splines_4.4.3              
##  [57] zlibbioc_1.52.0             parallel_4.4.3             
##  [59] restfulr_0.0.15             XVector_0.46.0             
##  [61] matrixStats_1.5.0           vctrs_0.6.5                
##  [63] Matrix_1.7-3                jsonlite_1.9.1             
##  [65] IRanges_2.40.1              S4Vectors_0.44.0           
##  [67] dendextend_1.19.0           adjclust_0.6.10            
##  [69] locfit_1.5-9.12             limma_3.62.2               
##  [71] jquerylib_0.1.4             glue_1.8.0                 
##  [73] codetools_0.2-20            stringi_1.8.4              
##  [75] gtable_0.3.6                GenomeInfoDb_1.42.3        
##  [77] GenomicRanges_1.58.0        BiocIO_1.16.0              
##  [79] UCSC.utils_1.2.0            munsell_0.5.1              
##  [81] tibble_3.2.1                pillar_1.10.1              
##  [83] rhdf5filters_1.18.0         csaw_1.40.0                
##  [85] htmltools_0.5.8.1           GenomeInfoDbData_1.2.13    
##  [87] BSgenome_1.74.0             R6_2.6.1                   
##  [89] sparseMatrixStats_1.18.0    diffHic_1.38.0             
##  [91] evaluate_1.0.3              lattice_0.22-6             
##  [93] Biobase_2.66.0              png_0.1-8                  
##  [95] Rsamtools_2.22.0            bslib_0.9.0                
##  [97] Rcpp_1.0.14                 InteractionSet_1.34.0      
##  [99] gridExtra_2.3               SparseArray_1.6.0          
## [101] xfun_0.51                   MatrixGenerics_1.18.0      
## [103] pkgconfig_2.0.3
```
