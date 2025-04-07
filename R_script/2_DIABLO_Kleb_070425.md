DIABLO - Identification of Signature Molecules of *K. pneumoniae* Serum
Responses
================

- [Load R packages](#load-r-packages)
- [Import all omics](#import-all-omics)
- [Data for MixOmics (DIABLO)](#data-for-mixomics-diablo)
- [Initial DIABLO model](#initial-diablo-model)
- [Tuning parameters](#tuning-parameters)
- [Tuning the number of features](#tuning-the-number-of-features)
- [Final DIABLO model](#final-diablo-model)
- [AUC for final model](#auc-for-final-model)
- [Extract data from final DIABLO
  model](#extract-data-from-final-diablo-model)
- [Data visualization](#data-visualization)
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

    library(igraph)

    ## 
    ## Attaching package: 'igraph'
    ## 
    ## The following objects are masked from 'package:lubridate':
    ## 
    ##     %--%, union
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union
    ## 
    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing
    ## 
    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     union

    set.seed(123) # for reproducibility

## Import all omics

    Kleb_trans_data <- read.csv("/Users/wmujchariyak/Desktop/Kleb_transcriptome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Kleb_prot_data <- read.csv("/Users/wmujchariyak/Desktop/Kleb_proteome_data.csv", header = TRUE) %>% column_to_rownames(var = "X")
    Kleb_met_GC_data <- read.csv("/Users/wmujchariyak/Desktop/Kleb_metabolome_GC_data.csv", header = TRUE) %>% column_to_rownames(var = "X")

## Data for MixOmics (DIABLO)

    data <- list(Transcript = Kleb_trans_data, 
                 Protein = Kleb_prot_data, 
                 Metabolite_GC = Kleb_met_GC_data)
    lapply(data, dim) # check their dimensions (24 rows)

    ## $Transcript
    ## [1]   24 5580
    ## 
    ## $Protein
    ## [1]   24 2013
    ## 
    ## $Metabolite_GC
    ## [1]  24 139

    Y <- factor(rep(c("RPMI", "SERUM"), 12)) # set the response variable as the Y df
    summary(Y)

    ##  RPMI SERUM 
    ##    12    12

## Initial DIABLO model

    design = matrix(0.1, ncol = length(data), nrow = length(data), 
                    dimnames = list(names(data), names(data)))
    diag(design) = 0 # set diagonal to 0s
    design

    ##               Transcript Protein Metabolite_GC
    ## Transcript           0.0     0.1           0.1
    ## Protein              0.1     0.0           0.1
    ## Metabolite_GC        0.1     0.1           0.0

## Tuning parameters

\*\* [LOO-CV
issue](https://mixomics-users.discourse.group/t/loo-cv-outputs-na-in-choice-ncomp/145/3) -
you can look at the perf plot and decide on the number of components
(Here, ncomp = 1).

    # Form basic DIABLO model with  an arbitrarily high number of components (ncomp = 5)
    basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design)

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

    # Tuning the number of components ("loo" = leave-one-out cross-validation, for small sample sizes)
    perf.diablo = perf(basic.diablo.model, validation = 'loo', nrepeat = 1) 
    plot(perf.diablo) # plot output of tuning

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_1.png)

    # Set the optimal ncomp value
    ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
    # Show the optimal choice for ncomp for each dist metric
    perf.diablo$choice.ncomp$WeightedVote 

    ## NULL

    # Check output - basic DIABLO model
    plotIndiv(basic.diablo.model, ind.names=FALSE, legend=TRUE)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_2.png)

    plotLoadings(basic.diablo.model)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_3.png)

## Tuning the number of features

\*\* This may take some time to run (~ 5 min)

    # Set grid of values for each component to test
    test.keepX = list (Transcript = seq(20,100,20), 
                       Protein = seq(20,100,20),
                       Metabolite_GC = seq(10,100,10)) 

    # Note: seq(10,100,10) = start from 10 and increase by 10 until 100

    # Run the feature selection tuning (ncomp = 1 - from tuning the number of components)
    tune.SA = tune.block.splsda(X = data, Y = Y, ncomp = 1, 
                                  test.keepX = test.keepX, design = design,
                                  validation = 'Mfold', folds = 5, nrepeat = 10,
                                  dist = "centroids.dist")

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

    ## 
    ## You have provided a sequence of keepX of length:  5 for block Transcript and  5 for block Protein and 10 for block Metabolite_GC.
    ## This results in 250 models being fitted for each component and each nrepeat, this may take some time to run, be patient!

    ## 
    ## You can look into the 'BPPARAM' argument to speed up computation time.

    # Set the optimal values of features to retain
    list.keepX = tune.SA$choice.keepX 
    list.keepX

    ## $Transcript
    ## [1] 20
    ## 
    ## $Protein
    ## [1] 20
    ## 
    ## $Metabolite_GC
    ## [1] 10

## Final DIABLO model

    # Final model (use arbitrary number for # variables from our preliminary results in "Tuning the number of features")

    list.keepX <-  list (Transcript =  20, 
                        Protein =  20,
                        Metabolite_GC = 10) 

    # Set the optimised DIABLO model
    final.diablo.model = block.splsda(X = data, Y = Y, ncomp = 1, 
                                      keepX = list.keepX, design = design)

    ## Design matrix has changed to include Y; each block will be
    ##             linked to Y.

    final.diablo.model$design # design matrix for the final model

    ##               Transcript Protein Metabolite_GC Y
    ## Transcript           0.0     0.1           0.1 1
    ## Protein              0.1     0.0           0.1 1
    ## Metabolite_GC        0.1     0.1           0.0 1
    ## Y                    1.0     1.0           1.0 0

## AUC for final model

    auc.diablo.trans.com <- auroc(final.diablo.model, roc.block = "Transcript", roc.comp = 1,
                            print = FALSE)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_4.png)

    auc.diablo.prot.com <- auroc(final.diablo.model, roc.block = "Protein", roc.comp = 1,
                            print = FALSE)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_5.png)

    auc.diablo.met.GC.com <- auroc(final.diablo.model, roc.block = "Metabolite_GC", roc.comp = ,
                            print = FALSE)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_6.png)

## Extract data from final DIABLO model

    # Variables selected on component 1
    trans_var <- selectVar(final.diablo.model, block = 'Transcript', comp = 1)
    prot_var <- selectVar(final.diablo.model, block = 'Protein', comp = 1)
    met_GC_var <- selectVar(final.diablo.model, block = 'Metabolite_GC', comp = 1)

    # Correlation matrix from the circos plot
    corMat <- circosPlot(final.diablo.model, cutoff = 0.7)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_7.png)

    corMat2 <- as.data.frame(corMat)

## Data visualization

    # Sample plots
    plotDiablo(final.diablo.model, ncomp = 1)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_8.png)

    # Plot Loadings
    plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median', size.name = 0.65)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_9.png)

    # Circos plot 
    circosPlot(final.diablo.model, cutoff = 0.9, comp = 1, line = TRUE, 
               color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
               color.cor = c("chocolate3","grey20"), size.labels = 1.2, size.variables = 0.5, size.legend = 1)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_10.png)

    circosPlot(final.diablo.model, cutoff = 0.95, comp = 1, line = TRUE, 
               color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
               color.cor = c("chocolate3","grey20"), size.labels = 1.2, size.variables = 0.5, size.legend = 1)

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_11.png)

## Save

    # Save data for further visualization in Cytoscape
    myNetwork <- network(final.diablo.model, blocks = c(1,2,3), cutoff = 0.9) 

![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_12.png)

    #write_graph(myNetwork$gR, file = "/Users/wmujchariyak/Desktop/myNetwork_conserved_Kleb.gml", format = "gml")

    myNetwork2 <- network(final.diablo.model, blocks = c(1,2,3), cutoff = 0.95) 
![](https://github.com/warasinee/Multiomics_Case_Study/blob/main/Image/Report_plots/2_DIABLO_plot_13.png)
    #write_graph(myNetwork2$gR, file = "/Users/wmujchariyak/Desktop/myNetwork_conserved_Kleb_cor095.gml", format = "gml")

