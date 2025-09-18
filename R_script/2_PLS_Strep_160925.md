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
    ## ✔ lubridate 1.9.4     ✔ tibble    3.3.0
    ## ✔ purrr     1.0.4     ✔ tidyr     1.3.1

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
                 Metabolite_GC = Strep_met_GC_data)
    lapply(data, dim) # check their dimensions (24 rows)

    ## $Transcript
    ## [1]   60 2239
    ## 
    ## $Protein
    ## [1]  60 993
    ## 
    ## $Metabolite_GC
    ## [1]  60 129

    Y <- factor(rep(c("RPMI", "SERUM"), 30)) # set the response variable as the Y df (based on conditions)
    Y2 <- factor(rep(c("5448_1", "HKU419","PS003","PS006","SP444"), each = 12)) # set the response variable as the Y df  (based on strains)

# Pairwise PLS Comparisons (ncomp=1)

\*\* Calculate correlation between omics datasets which can be
beneficial for setting of design matrix

    # generate pairwise PLS models 
    pls1 <- pls(data$Transcript, data$Protein, ncomp = 1) 
    pls2 <- pls(data$Protein, data$Metabolite_GC, ncomp = 1)
    pls3 <- pls(data$Transcript, data$Metabolite_GC, ncomp = 1)

    # calculate correlation 
    cor(pls1$variates$X, pls1$variates$Y) #cor trans vs prot = 0.9319566 

    ##           comp1
    ## comp1 0.9319566

    cor(pls2$variates$X, pls2$variates$Y) #cor prot vs met_GC = 0.9127586

    ##           comp1
    ## comp1 0.9127586

    cor(pls3$variates$X, pls3$variates$Y) #cor trans vs met_GC = 0.9343444

    ##           comp1
    ## comp1 0.9343444

## sPLS-DA: based on conditions

\*\* An sPLS-DA analysis will help refine the sample clusters and select
a small subset of variables relevant to discriminate each class.

    ## Transcript 
    # Load the data
    X.trans <- data$Transcript
    dim(X.trans); length(Y)

    ## [1]   60 2239

    ## [1] 60

    # Initial model
    plsda.trans <- plsda(X.trans,Y, ncomp = 10)

    # Tuning number of components
    perf.plsda.trans <- perf(plsda.trans, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.trans$error.rate) 

    ## $overall
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.01666667     0.01666667       0.01666667
    ## comp2  0.00000000     0.00000000       0.00000000
    ## comp3  0.00000000     0.00000000       0.00000000
    ## comp4  0.00000000     0.00000000       0.00000000
    ## comp5  0.00000000     0.00000000       0.00000000
    ## comp6  0.00000000     0.00000000       0.00000000
    ## comp7  0.00000000     0.00000000       0.00000000
    ## comp8  0.00000000     0.00000000       0.00000000
    ## comp9  0.00000000     0.00000000       0.00000000
    ## comp10 0.00000000     0.00000000       0.00000000
    ## 
    ## $BER
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.01666667     0.01666667       0.01666667
    ## comp2  0.00000000     0.00000000       0.00000000
    ## comp3  0.00000000     0.00000000       0.00000000
    ## comp4  0.00000000     0.00000000       0.00000000
    ## comp5  0.00000000     0.00000000       0.00000000
    ## comp6  0.00000000     0.00000000       0.00000000
    ## comp7  0.00000000     0.00000000       0.00000000
    ## comp8  0.00000000     0.00000000       0.00000000
    ## comp9  0.00000000     0.00000000       0.00000000
    ## comp10 0.00000000     0.00000000       0.00000000

    plot(perf.plsda.trans, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-1.png)

    # The error rate decreases and reaches a minimum for ncomp = 2 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.trans <- plsda(X.trans,Y, ncomp = 2)

    plotIndiv_trans <- plotIndiv(final.plsda.trans, ind.names = T, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Transcriptome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-2.png)

    ######################################
    ## Proteomics 
    # Load the data
    X.prot <- data$Protein
    dim(X.prot); length(Y)

    ## [1]  60 993

    ## [1] 60

    # Initial model
    plsda.prot <- plsda(X.prot,Y, ncomp = 10)

    # Tuning number of components
    perf.plsda.prot <- perf(plsda.prot, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.prot$error.rate) 

    ## $overall
    ##        max.dist centroids.dist mahalanobis.dist
    ## comp1         0              0                0
    ## comp2         0              0                0
    ## comp3         0              0                0
    ## comp4         0              0                0
    ## comp5         0              0                0
    ## comp6         0              0                0
    ## comp7         0              0                0
    ## comp8         0              0                0
    ## comp9         0              0                0
    ## comp10        0              0                0
    ## 
    ## $BER
    ##        max.dist centroids.dist mahalanobis.dist
    ## comp1         0              0                0
    ## comp2         0              0                0
    ## comp3         0              0                0
    ## comp4         0              0                0
    ## comp5         0              0                0
    ## comp6         0              0                0
    ## comp7         0              0                0
    ## comp8         0              0                0
    ## comp9         0              0                0
    ## comp10        0              0                0

    plot(perf.plsda.prot, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-3.png)

    # The error rate decreases and reaches a minimum for ncomp = 1 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.prot <- plsda(X.prot,Y, ncomp = 2) # we use ncomp =2 so that we can visualize with plotIndiv

    plotIndiv_prot <- plotIndiv(final.plsda.prot, ind.names = T, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Proteome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-4.png)

    ######################################
    ## Met - sPLS-DA
    # Load the data
    X.met <- data$Metabolite_GC
    dim(X.met); length(Y)

    ## [1]  60 129

    ## [1] 60

    # Initial model
    plsda.met <- plsda(X.met,Y, ncomp = 10)

    # Tuning number of components
    perf.plsda.met <- perf(plsda.met, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.met$error.rate) 

    ## $overall
    ##        max.dist centroids.dist mahalanobis.dist
    ## comp1         0              0                0
    ## comp2         0              0                0
    ## comp3         0              0                0
    ## comp4         0              0                0
    ## comp5         0              0                0
    ## comp6         0              0                0
    ## comp7         0              0                0
    ## comp8         0              0                0
    ## comp9         0              0                0
    ## comp10        0              0                0
    ## 
    ## $BER
    ##        max.dist centroids.dist mahalanobis.dist
    ## comp1         0              0                0
    ## comp2         0              0                0
    ## comp3         0              0                0
    ## comp4         0              0                0
    ## comp5         0              0                0
    ## comp6         0              0                0
    ## comp7         0              0                0
    ## comp8         0              0                0
    ## comp9         0              0                0
    ## comp10        0              0                0

    plot(perf.plsda.met, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-5.png)

    # The error rate decreases and reaches a minimum for ncomp = 1 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.met <- plsda(X.met,Y, ncomp = 2) # we use ncomp =2 so that we can visualize with plotIndiv

    plotIndiv_met <- plotIndiv(final.plsda.met, ind.names = T, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Metabolome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-5-6.png)

## sPLS-DA: based on strains

    ## Transcript 
    # Load the data
    X.trans <- data$Transcript
    dim(X.trans); length(Y2)

    ## [1]   60 2239

    ## [1] 60

    # Initial model
    plsda.trans.2 <- plsda(X.trans,Y2, ncomp = 10)

    # Tuning number of components
    perf.plsda.trans.2 <- perf(plsda.trans.2, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.trans.2$error.rate) 

    ## $overall
    ##         max.dist centroids.dist mahalanobis.dist
    ## comp1  0.6000000      0.1833333       0.18333333
    ## comp2  0.6166667      0.0500000       0.05000000
    ## comp3  0.3500000      0.0000000       0.01666667
    ## comp4  0.0000000      0.0000000       0.00000000
    ## comp5  0.0000000      0.0000000       0.00000000
    ## comp6  0.0000000      0.0000000       0.00000000
    ## comp7  0.0000000      0.0000000       0.00000000
    ## comp8  0.0000000      0.0000000       0.00000000
    ## comp9  0.0000000      0.0000000       0.00000000
    ## comp10 0.0000000      0.0000000       0.00000000
    ## 
    ## $BER
    ##         max.dist centroids.dist mahalanobis.dist
    ## comp1  0.6000000      0.1833333       0.18333333
    ## comp2  0.6166667      0.0500000       0.05000000
    ## comp3  0.3500000      0.0000000       0.01666667
    ## comp4  0.0000000      0.0000000       0.00000000
    ## comp5  0.0000000      0.0000000       0.00000000
    ## comp6  0.0000000      0.0000000       0.00000000
    ## comp7  0.0000000      0.0000000       0.00000000
    ## comp8  0.0000000      0.0000000       0.00000000
    ## comp9  0.0000000      0.0000000       0.00000000
    ## comp10 0.0000000      0.0000000       0.00000000

    plot(perf.plsda.trans.2, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-1.png)

    # The error rate decreases and reaches a minimum for ncomp = 4 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.trans.2 <- plsda(X.trans,Y2, ncomp = 4)

    plotIndiv_trans.2 <- plotIndiv(final.plsda.trans.2, ind.names = F, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Transcriptome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-2.png)

    ######################################
    ## Proteomics 
    # Load the data
    X.prot <- data$Protein
    dim(X.prot); length(Y2)

    ## [1]  60 993

    ## [1] 60

    # Initial model
    plsda.prot.2 <- plsda(X.prot,Y2, ncomp = 10)

    # Tuning number of components
    perf.plsda.prot.2 <- perf(plsda.prot.2, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.prot.2$error.rate) 

    ## $overall
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.60000000     0.25000000       0.25000000
    ## comp2  0.40000000     0.11666667       0.03333333
    ## comp3  0.20000000     0.11666667       0.05000000
    ## comp4  0.03333333     0.03333333       0.00000000
    ## comp5  0.00000000     0.00000000       0.00000000
    ## comp6  0.00000000     0.00000000       0.00000000
    ## comp7  0.00000000     0.00000000       0.00000000
    ## comp8  0.00000000     0.00000000       0.00000000
    ## comp9  0.00000000     0.00000000       0.00000000
    ## comp10 0.00000000     0.00000000       0.00000000
    ## 
    ## $BER
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.60000000     0.25000000       0.25000000
    ## comp2  0.40000000     0.11666667       0.03333333
    ## comp3  0.20000000     0.11666667       0.05000000
    ## comp4  0.03333333     0.03333333       0.00000000
    ## comp5  0.00000000     0.00000000       0.00000000
    ## comp6  0.00000000     0.00000000       0.00000000
    ## comp7  0.00000000     0.00000000       0.00000000
    ## comp8  0.00000000     0.00000000       0.00000000
    ## comp9  0.00000000     0.00000000       0.00000000
    ## comp10 0.00000000     0.00000000       0.00000000

    plot(perf.plsda.prot.2, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-3.png)

    # The error rate decreases and reaches a minimum for ncomp = 5 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.prot.2 <- plsda(X.prot,Y2, ncomp = 5) 

    plotIndiv_prot.2 <- plotIndiv(final.plsda.prot.2, ind.names = F, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Proteome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-4.png)

    ######################################
    ## Met - sPLS-DA
    # Load the data
    X.met <- data$Metabolite_GC
    dim(X.met); length(Y2)

    ## [1]  60 129

    ## [1] 60

    # Initial model
    plsda.met.2 <- plsda(X.met,Y2, ncomp = 10)

    # Tuning number of components
    perf.plsda.met.2 <- perf(plsda.met.2, validation = 'loo', 
                      progressBar = FALSE,  # Set to TRUE to track progress
                      nrepeat = 1)    

    head(perf.plsda.met.2$error.rate) 

    ## $overall
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.80000000     0.70000000       0.70000000
    ## comp2  0.45000000     0.31666667       0.31666667
    ## comp3  0.26666667     0.10000000       0.15000000
    ## comp4  0.21666667     0.21666667       0.08333333
    ## comp5  0.08333333     0.06666667       0.05000000
    ## comp6  0.05000000     0.05000000       0.03333333
    ## comp7  0.03333333     0.05000000       0.03333333
    ## comp8  0.03333333     0.01666667       0.01666667
    ## comp9  0.01666667     0.01666667       0.03333333
    ## comp10 0.01666667     0.01666667       0.00000000
    ## 
    ## $BER
    ##          max.dist centroids.dist mahalanobis.dist
    ## comp1  0.80000000     0.70000000       0.70000000
    ## comp2  0.45000000     0.31666667       0.31666667
    ## comp3  0.26666667     0.10000000       0.15000000
    ## comp4  0.21666667     0.21666667       0.08333333
    ## comp5  0.08333333     0.06666667       0.05000000
    ## comp6  0.05000000     0.05000000       0.03333333
    ## comp7  0.03333333     0.05000000       0.03333333
    ## comp8  0.03333333     0.01666667       0.01666667
    ## comp9  0.01666667     0.01666667       0.03333333
    ## comp10 0.01666667     0.01666667       0.00000000

    plot(perf.plsda.met.2, sd = TRUE, legend.position = 'horizontal')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-5.png)

    # The error rate decreases and reaches a minimum for ncomp = 5 for the max.dist distance. These parameters will be included in further analyses.

    # Final model
    final.plsda.met.2 <- plsda(X.met,Y2, ncomp = 5) # we use ncomp =2 so that we can visualize with plotIndiv

    plotIndiv_met.2 <- plotIndiv(final.plsda.met.2, ind.names = F, legend=TRUE,
              comp=c(1,2), ellipse = TRUE, 
              title = 'PLS-DA on Metabolome comp 1-2',
              X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

![](2_PLS_Strep_160925_files/figure-markdown_strict/unnamed-chunk-6-6.png)
\## Save

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_trans.pdf", height = 5, width = 8)
    #plotIndiv_trans
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_prot.pdf", height = 5, width = 8)
    #plotIndiv_prot
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_met_GC.pdf", height = 5, width = 8)
    #plotIndiv_met
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_trans2.pdf", height = 5, width = 8)
    #plotIndiv_trans.2
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_prot2.pdf", height = 5, width = 8)
    #plotIndiv_prot.2
    #dev.off()

    #pdf("/Users/wmujchariyak/Desktop/Strep_plotIndiv_plsda_met_GC2.pdf", height = 5, width = 8)
    #plotIndiv_met.2
    #dev.off()
