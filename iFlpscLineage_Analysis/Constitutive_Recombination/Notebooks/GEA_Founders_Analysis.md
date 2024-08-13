GEA_Founders_Analysis
================
Alvaro Regano
27/07/2023

# GEA Founders Analysis

# Objective

Take a look at the new GEA line and conclude that it is an improvement
on the last by checking for higher expression of iFlpscLineage BCs in
the new founder as well as a higher number of cells expressing the
transgene

# Analysis

The dataset conditions is as follows

    ##   Founder Recomb Promoter Origin CMO
    ## 1   GEA12   FULL    Super   GEA3 301
    ## 2   GEA12   NONE    Super   GEA3 302
    ## 3   GEA18   FULL      R26   GEA3 303
    ## 4    GEA3   NONE      R26   GEA3 304
    ## 5   GEA27   NONE      R26   GEA2 305

    ## 
    ## GEA12_FULL_recomb_iFlpscL   GEA12_NO_recomb_iFlpscL GEA18_FULL_recomb_iFlpscL    GEA3_NO_recomb_iFlpscL   GEA27_NO_recomb_iFlpscL 
    ##                       802                       912                       556                       552                       534

## Module Scores of the iFlpscLineage cassette

First We will look at the Barcodes grouped to see how the expression
rates change

    ## png 
    ##   2

    ## png 
    ##   2

Violin Plots show that GEA12 with Full recombination has a remarkable
level of expression of all iFlpscLineage Arrays

<img src="GEA_Founders_Analysis_files/figure-gfm/UMAPwoInt3-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/UMAPwoInt3-2.png" width="200%" style="display: block; margin: auto;" />

    ## TableGrob (5 x 5) "guide-box": 2 grobs
    ##   z     cells                  name           grob
    ## 1 1 (3-3,3-3)                guides gtable[layout]
    ## 2 0 (2-4,2-4) legend.box.background zeroGrob[NULL]

<img src="GEA_Founders_Analysis_files/figure-gfm/UMAPwoInt3-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/UMAPwoInt3-4.png" width="200%" style="display: block; margin: auto;" />

    ## png 
    ##   2

    ## png 
    ##   2

We can also look at the % of cells expressing iFlpscLineage in the
dataset

## Looking at Expression rates of iFlpscLinage per cell

<img src="GEA_Founders_Analysis_files/figure-gfm/UMAPCondSeuratInt-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/UMAPCondSeuratInt-2.png" width="200%" style="display: block; margin: auto;" />

Now in terms of % of counts in the dataset

<img src="GEA_Founders_Analysis_files/figure-gfm/barplot counts-1.png" width="200%" style="display: block; margin: auto;" />

Surprisingly there are some cells in the GEA12_FULL_recomb that seem to
have a huge rate of expression of FRT_BCs (over 60%)

We can take a look at rates of expression among arrays

<img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per array barplots-1.png" width="200%" style="display: block; margin: auto;" />

Let’s check the number and percentage of cells expressing more than one
Array in the dataset

<img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per array-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per array-2.png" width="200%" style="display: block; margin: auto;" />

For the GEA12 Founder the proportion of cells expressing iFlpscLineage
Barcodes is huge, over 80% of them express at least one iFlpscLineage
Array.

## Looking at the Distribution of FRT-BCs throughout the populations

### Looking at FRT-BC Mean Counts per cell

<img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per condition-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per condition-2.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per condition-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs per condition-4.png" width="200%" style="display: block; margin: auto;" />

### Looking at FRT-BC Counts

And Now instead of Mean counts per cell lets look at the total counts in
the dataset

<img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs counts per condition-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs counts per condition-2.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs counts per condition-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/FRT_BCs counts per condition-4.png" width="200%" style="display: block; margin: auto;" />

## Doing a ratio analysis between Fwd and Inverted BC segment

<img src="GEA_Founders_Analysis_files/figure-gfm/BC ratios per condition plot-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/BC ratios per condition plot-2.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/BC ratios per condition plot-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/BC ratios per condition plot-4.png" width="200%" style="display: block; margin: auto;" />
The ratio analysis between Inv/Fwd BC segments indicate a good level of
Iverted BCs present in the dataset

## Looking at the Distribution of other transgenes in the populations

<img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes per condition-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes per condition-2.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes per condition-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes per condition-4.png" width="200%" style="display: block; margin: auto;" />

### Looking at Transgene Counts

And Now instead of Mean counts per cell lets look at the total counts in
the dataset

<img src="GEA_Founders_Analysis_files/figure-gfm/Transgene counts per condition-1.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgene counts per condition-2.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgene counts per condition-3.png" width="200%" style="display: block; margin: auto;" /><img src="GEA_Founders_Analysis_files/figure-gfm/Transgene counts per condition-4.png" width="200%" style="display: block; margin: auto;" />

## Looking at the other Transgenes present in the dataset

<img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes-1.png" width="200%" style="display: block; margin: auto;" />

<img src="GEA_Founders_Analysis_files/figure-gfm/Transgenes VlnPlot-1.png" width="200%" style="display: block; margin: auto;" />

# Conclusion

GEA12 undergoing full recombination is by far the best of all conditions
present. Induction has worked by seeing the levels of expression of FlpE
are only high in Fully recombined conditions There are issues in the
correct alignment of Barcodes as in non recombined conditions appear
levels of Reverse Barcodes. However, we can get rid of those issues by
doing Nanopore sequencing, as Irepan’s RFR_ONT pipeline can correctly
address that only properly oriented full transcript iFlpscLineage
Barcodes pass the Analysis. Here I can only select for the most frequent
Barcode in each cell, but I have not done it as it can bias the counts,
favoring cells with degenerate BCs even more, and overall analysis.

# R session

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 24.04 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Madrid
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] cowplot_1.1.3      RColorBrewer_1.1-3 tidyr_1.3.1        patchwork_1.2.0    yaml_2.3.10        rmarkdown_2.27     dplyr_1.1.4        ggplot2_3.5.1      Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4           knitr_1.48        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.4            magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.3-2    matrixStats_1.3.0      ggridges_0.5.6         compiler_4.3.3         png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.2.0          labeling_0.4.3         utf8_1.2.4             promises_1.3.0         purrr_1.0.2            xfun_0.46              jsonlite_1.8.8         goftest_1.2-3          highr_0.11             later_1.3.2            spatstat.utils_3.0-5   irlba_2.3.5.1          parallel_4.3.3         cluster_2.1.6          R6_2.5.1               ica_1.0-3              stringi_1.8.4          spatstat.data_3.1-2    reticulate_1.38.0      parallelly_1.38.0      spatstat.univar_3.0-0  lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.13            tensor_1.5             future.apply_1.11.2    zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.15          Matrix_1.6-5           splines_4.3.3          igraph_2.0.3           tidyselect_1.2.1       abind_1.4-5            rstudioapi_0.16.0      spatstat.random_3.3-1  codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.3-1 listenv_0.9.1          lattice_0.22-5         tibble_3.2.1           plyr_1.8.9             withr_3.0.1            shiny_1.9.1            ROCR_1.0-11            evaluate_0.24.0        Rtsne_0.17             future_1.34.0          fastDummies_1.7.3      survival_3.5-8         polyclip_1.10-7        fitdistrplus_1.2-1     pillar_1.9.0           KernSmooth_2.23-22     plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         munsell_0.5.1          scales_1.3.0           globals_0.16.3         xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2         tools_4.3.3            data.table_1.15.4      RSpectra_0.16-2        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        grid_4.3.3             colorspace_2.1-1       nlme_3.1-164           cli_3.6.3              spatstat.sparse_3.1-0  spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2      uwot_0.2.2             gtable_0.3.5           digest_0.6.36          progressr_0.14.0       ggrepel_0.9.5          farver_2.1.2           htmlwidgets_1.6.4      htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.12              MASS_7.3-60.0.1
