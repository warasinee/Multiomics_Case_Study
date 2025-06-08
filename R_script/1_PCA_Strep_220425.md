Unsupervised Analysis - PCA
================

- [Load R packages](#load-r-packages)
- [Import all omics](#import-all-omics)
- [Data for MixOmics](#data-for-mixomics)
- [PCA](#pca)
- [Data visualization - Preliminary](#data-visualization---preliminary)
- [Data visualization - Final](#data-visualization---final)
- [Save](#save)


## Load R packages
  

library(mixOmics)


    ## Loading required package: MASS

    ## Loading required package: lattice

    ## Loading required package: ggplot2

    ## 
    ## Loaded mixOmics 6.22.0
    ## Thank you for using mixOmics!
    ## Tutorials: http://mixomics.org
    ## Bookdown vignette: https://mixomicsteam.github.io/Bookdown
    ## Questions, issues: Follow the prompts at http://mixomics.org/contact-us
    ## Cite us:  citation('mixOmics')


library(tidyverse)


    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ purrr::map()    masks mixOmics::map()
    ## ✖ dplyr::select() masks MASS::select()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors


set.seed(123) # for reproducibility


## Import all omics


    Strep_trans_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_transcriptome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Strep_prot_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_proteome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Strep_met_GC_data <- read.csv("/Users/wmujchariyak/Desktop/Strep_metabolome_GC_data.csv", header = TRUE) %>% column_to_rownames(var = "X")


## Data for MixOmics


    data <- list(Transcript = Strep_trans_data, 
                 Protein = Strep_prot_data, 
                 Metabolite_LC = Strep_met_GC_data)
    lapply(data, dim) # check their dimensions (24 rows)

    ## $Transcript
    ## [1]   60 2239
    ## 
    ## $Protein
    ## [1]  60 993
    ## 
    ## $Metabolite_LC
    ## [1]  60 129

    Y <- factor(rep(c("RPMI", "SERUM"), 30)) # set the response variable as the Y df
    summary(Y)

    ##  RPMI SERUM 
    ##    30    30

## PCA

    data_name <- c("Transcript", "Protein", "Metabolite_LC")
    list.final.pca <- list()
    plot.tune.pca <- list()
    pca.cum.var <- list()

    for (i in seq_along(data_name)){
      set.seed(123) # for reproducibility
      data_pca <- data[[i]]
      # Choosing number of components
      tune.pca <- tune.pca(data_pca, ncomp = 10, scale = TRUE)
      plot.tune.pca[[i]] <- tune.pca
      # Outputs cumulative proportion of variance
      pca.cum.var[[i]] <- tune.pca$cum.var       
      # Final pca 
      final.pca <- pca(data_pca, ncomp = 3, scale = TRUE)
      list.final.pca[[i]] <- final.pca
    }

    # Note: Results from "tune.pca" suggest that using 3 components can cover almost 50% of variation in the datasets

## Data visualization - Preliminary

    # Trans
    plotLoadings(list.final.pca[[1]], comp = 1, method = 'mean', contrib = 'max', size.title = rel(1))

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/1_PCA_plot_1.png)

    trans.pca.var <- selectVar(list.final.pca[[1]], comp = 1)$value
    trans.pca.var$Rank <- rank(abs(trans.pca.var$value.var))

    # Prot
    plotLoadings(list.final.pca[[2]], comp = 1, method = 'mean', contrib = 'max', size.title = rel(1))

![](1_PCA_Strep_220425_files/figure-markdown_strict/unnamed-chunk-5-2.png)

    prot.ms1.pca.var <- selectVar(list.final.pca[[2]], comp = 1)$value
    prot.ms1.pca.var$Rank <- rank(abs(prot.ms1.pca.var$value.var))

    # Met - LC
    plotLoadings(list.final.pca[[3]], comp = 1, method = 'mean', contrib = 'max', size.title = rel(1))

![](1_PCA_Strep_220425_files/figure-markdown_strict/unnamed-chunk-5-3.png)

    met.LC.pca.var <- selectVar(list.final.pca[[3]], comp = 1)$value
    met.LC.pca.var$Rank <- rank(abs(met.LC.pca.var$value.var))

## Data visualization - Final

    plotIndiv_pca_trans <- plotIndiv(list.final.pca[[1]], group = Y, ind.names = FALSE,
                                      legend = TRUE, ellipse = T, style = "ggplot2", size.title = rel(2.5),
                                      size.xlabel = rel(1.5), size.ylabel = rel(1.5),
                                      title = 'Transcriptome') 

![](1_PCA_Strep_220425_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    plotIndiv_pca_prot <- plotIndiv(list.final.pca[[2]], group = Y, ind.names = FALSE,
                                     legend = TRUE, ellipse = T, style = "ggplot2", 
                                     size.title = rel(2.5),
                                     size.xlabel = rel(1.5), size.ylabel = rel(1.5),
                                     title = 'Proteome')

![](1_PCA_Strep_220425_files/figure-markdown_strict/unnamed-chunk-6-2.png)

    plotIndiv_pca_met_LC <- plotIndiv(list.final.pca[[3]], group = Y, ind.names = FALSE,
                                       legend = TRUE, ellipse = T, style = "ggplot2", 
                                       size.title = rel(2.5),
                                       size.xlabel = rel(1.5), size.ylabel = rel(1.5),
                                       title = 'Metabolome (LC-MS)')

![](1_PCA_Strep_220425_files/figure-markdown_strict/unnamed-chunk-6-3.png)

## Save

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_pca_trans.pdf", height = 5, width = 8)
    #plotIndiv_pca_trans
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_pca_prot.pdf", height = 5, width = 8)
    #plotIndiv_pca_prot
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_pca_met_LC.pdf", height = 5, width = 8)
    #plotIndiv_pca_met_LC
    #dev.off()
