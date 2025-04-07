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
