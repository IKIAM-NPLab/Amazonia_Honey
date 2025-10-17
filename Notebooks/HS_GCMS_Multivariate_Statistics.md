Chemical characterization of volatile compounds and biological activity
of stingless bee honeys from the Ecuadorian Amazon using Headspace Gas
Chromatography-Mass Spectrometry - Statistical analysis
================
Jefferson Pastuña
2024-04-17

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#before-to-start" id="toc-before-to-start">Before to start</a>
- <a href="#notame-workflow" id="toc-notame-workflow">Notame workflow</a>
  - <a href="#preprocessing" id="toc-preprocessing">Preprocessing</a>
  - <a href="#processing" id="toc-processing">Processing</a>
- <a href="#principal-component-analysis-pca"
  id="toc-principal-component-analysis-pca">Principal component analysis
  (PCA)</a>
  - <a href="#score-pca" id="toc-score-pca">Score PCA</a>
  - <a href="#loading-pca" id="toc-loading-pca">Loading PCA</a>
- <a href="#permanova-analysis" id="toc-permanova-analysis">PERMANOVA
  analysis</a>
  - <a href="#t-angustula" id="toc-t-angustula"><em>T. angustula</em></a>
  - <a href="#m-fasciculata" id="toc-m-fasciculata"><em>M.
    fasciculata</em></a>
  - <a href="#m-fuscopilosa" id="toc-m-fuscopilosa"><em>M.
    fuscopilosa</em></a>
- <a href="#heatmap-with-hca" id="toc-heatmap-with-hca">Heatmap with
  HCA</a>
- <a href="#antimicrobial-activity-correlation"
  id="toc-antimicrobial-activity-correlation">Antimicrobial activity
  correlation</a>
  - <a href="#using-the-pca" id="toc-using-the-pca">Using the PCA</a>
    - <a href="#score-pca-1" id="toc-score-pca-1">Score PCA</a>
    - <a href="#loading-pca-1" id="toc-loading-pca-1">Loading PCA</a>
  - <a href="#using-the-pearson-correlation"
    id="toc-using-the-pearson-correlation">Using the Pearson correlation</a>
    - <a href="#normality-of-the-antimicrobial-data"
      id="toc-normality-of-the-antimicrobial-data">Normality of the
      antimicrobial data</a>
    - <a href="#metabolomics-data-transformation"
      id="toc-metabolomics-data-transformation">Metabolomics data
      transformation</a>
    - <a href="#correlation-of-identified-metabolites"
      id="toc-correlation-of-identified-metabolites">Correlation of identified
      metabolites</a>
    - <a href="#correlation-of-unknown-metabolites"
      id="toc-correlation-of-unknown-metabolites">Correlation of unknown
      metabolites</a>
- <a href="#upset-plot" id="toc-upset-plot">UpSet plot</a>

# Introduction

This R Script aims to record the procedure for analyzing the volatile
profile of honeys from three different stingless bee species (*Melipona
fasciculata*, *Melipona fuscopilosa*, and *Tetragonisca angustula*).
Each step has a brief explanation, as well as code and graphics.

The data preprocessing workflow used was taken from the research paper
[“notame”: Workflow for Non-Targeted LC–MS Metabolic
Profiling](https://doi.org/10.3390/metabo10040135), which offers a wide
variety of functions for performing metabolomic profile analysis.

# Before to start

The “*notame*” package accepts a feature table that can be obtained
through different software, such as MZmine, MS-DIAL, and others, as
input. The feature list table was extracted using the MZmine software in
this case. The exported (.csv) file from MZmine was fixed to get the
final feature table input according to the “*notame*” package guideline.

The raw (.csv) file modifications can be summarized by adding and
renaming columns. The added columns “Column” and “Ion Mode” allow for
the analysis of samples with different types of columns and different
ionization polarities, respectively. Also, the cells corresponding to
mass (*m/z*) and retention time (*rt*) must be renamed so the “*notame*”
package can detect and process them.

# Notame workflow

The “*notame*” package and other dependency packages were installed as a
first step for the analysis.

``` r
# Notame package installation
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("antonvsdata/notame", ref = "v0.3.1")

# Notame library call
library(notame)
```

    ## Warning: package 'BiocGenerics' was built under R version 4.5.0

``` r
# Dependency packages installation
install_dependencies
```

    ## function (preprocessing = TRUE, extra = FALSE, batch_corr = FALSE, 
    ##     misc = FALSE, ...) 
    ## {
    ##     core_cran <- c("BiocManager", "cowplot", "missForest", "openxlsx", 
    ##         "randomForest", "RColorBrewer", "Rtsne")
    ##     core_bioconductor <- "pcaMethods"
    ##     extra_cran <- c("car", "doParallel", "devEMF", "ggbeeswarm", 
    ##         "ggdendro", "ggrepel", "ggtext", "Hmisc", "hexbin", "igraph", 
    ##         "lme4", "lmerTest", "MuMIn", "PERMANOVA", "PK", "rmcorr")
    ##     extra_bioconductor <- c("mixOmics", "supraHex")
    ##     extra_gitlab <- "CarlBrunius/MUVR"
    ##     batch_cran <- "fpc"
    ##     batch_bioconductor <- "RUVSeq"
    ##     batch_github <- NULL
    ##     batch_gitlab <- "CarlBrunius/batchCorr"
    ##     misc_cran <- c("knitr", "rmarkdown", "testthat")
    ##     if (preprocessing) {
    ##         install_helper(cran = core_cran, bioconductor = core_bioconductor, 
    ##             ...)
    ##     }
    ##     if (extra) {
    ##         install_helper(cran = extra_cran, bioconductor = extra_bioconductor, 
    ##             gitlab = extra_gitlab, ...)
    ##     }
    ##     if (batch_corr) {
    ##         install_helper(cran = batch_cran, bioconductor = batch_bioconductor, 
    ##             github = batch_github, gitlab = batch_gitlab, ...)
    ##     }
    ##     if (misc) {
    ##         install_helper(cran = misc_cran, ...)
    ##     }
    ## }
    ## <bytecode: 0x0000000010d50ef0>
    ## <environment: namespace:notame>

Then, a primary path and a log system were added to have a record of
each process executed.

``` r
# Main path
ppath <- "../Amazonia_Honey/"
# Log system
init_log(log_file = paste0(ppath, "../Result/notame_Result/HS_GCMS/HS_GCMS_log.txt"))
```

    ## INFO [2025-10-13 08:26:56] Starting logging

Next, the MZmine feature list table in “*notame*” format was loaded.

``` r
data <- read_from_excel(file = "../Data/Data_to_notame/HS_GCMS_Data_to_notame.xlsx",
                        sheet = 2, corner_row = 18, corner_column = "L",
                        split_by = c("Column", "Ion Mode"))
```

    ## INFO [2025-10-13 08:26:56] Corner detected correctly at row 18, column L
    ## INFO [2025-10-13 08:26:56] 
    ## Extracting sample information from rows 1 to 18 and columns M to AV
    ## INFO [2025-10-13 08:26:56] Replacing spaces in sample information column names with underscores (_)
    ## INFO [2025-10-13 08:26:56] Naming the last column of sample information "Datafile"
    ## INFO [2025-10-13 08:26:56] 
    ## Extracting feature information from rows 19 to 704 and columns A to L
    ## INFO [2025-10-13 08:26:56] Creating Split column from Column, Ion Mode
    ## INFO [2025-10-13 08:26:56] Feature_ID column not found, creating feature IDs
    ## INFO [2025-10-13 08:26:56] Identified m/z column mass and retention time column RT
    ## INFO [2025-10-13 08:26:56] Identified m/z column mass and retention time column RT
    ## INFO [2025-10-13 08:26:56] Creating feature IDs from Split, m/z and retention time
    ## INFO [2025-10-13 08:26:56] Replacing dots (.) in feature information column names with underscores (_)
    ## INFO [2025-10-13 08:26:56] 
    ## Extracting feature abundances from rows 19 to 704 and columns M to AV
    ## INFO [2025-10-13 08:26:56] 
    ## Checking sample information
    ## INFO [2025-10-13 08:26:56] QC column generated from rows containing 'QC'
    ## INFO [2025-10-13 08:26:56] Sample ID autogenerated from injection orders and prefix ID_
    ## INFO [2025-10-13 08:26:56] Checking that feature abundances only contain numeric values
    ## INFO [2025-10-13 08:26:56] 
    ## Checking feature information
    ## INFO [2025-10-13 08:26:56] Checking that feature IDs are unique and not stored as numbers
    ## INFO [2025-10-13 08:26:56] Checking that m/z and retention time values are reasonable
    ## INFO [2025-10-13 08:26:56] Identified m/z column mass and retention time column RT
    ## INFO [2025-10-13 08:26:56] Identified m/z column mass and retention time column RT

Once the data was loaded, the next step was to create a MetaboSet to
work with R objects from now on.

``` r
modes <- construct_metabosets(exprs = data$exprs, 
                              pheno_data = data$pheno_data, 
                              feature_data = data$feature_data,
                              group_col = "Group")
```

    ## Initializing the object(s) with unflagged features
    ## INFO [2025-10-13 08:26:56] 
    ## Checking feature information
    ## INFO [2025-10-13 08:26:57] Checking that feature IDs are unique and not stored as numbers
    ## INFO [2025-10-13 08:26:57] Checking that feature abundances only contain numeric values
    ## INFO [2025-10-13 08:26:57] Setting row and column names of exprs based on feature and pheno data

## Preprocessing

The first step of the preprocessing is to change the features with a
value equal to 0 to NA.

``` r
# Data extraction
mode <- modes$RTX5MS_EI
# Change 0 value to NA
mode <- mark_nas(mode, value = 0)
```

Then, features with low detection rates are flagged and can be ignored
or removed in subsequent analysis. The “*notame*” package uses two
criteria to flag these features: the feature’s presence in a percentage
of QC injections and the feature’s presence in a percentage within a
sample group or class.

``` r
# Low detection rate
mode <- flag_detection(mode, qc_limit = 7/9, group_limit = 2/3)
```

    ## INFO [2025-10-13 08:26:57] 
    ## 2% of features flagged for low detection rate

The following preprocessing step is drift correction, which is applied
using smoothed cubic spline regression.

``` r
# Drift correction
corrected <- correct_drift(mode)
```

    ## INFO [2025-10-13 08:26:57] 
    ## Starting drift correction at 2025-10-13 08:26:57.662858
    ## INFO [2025-10-13 08:26:59] Drift correction performed at 2025-10-13 08:26:59.333954
    ## INFO [2025-10-13 08:27:00] Inspecting drift correction results 2025-10-13 08:27:00.198003
    ## INFO [2025-10-13 08:27:01] Drift correction results inspected at 2025-10-13 08:27:01.853098
    ## INFO [2025-10-13 08:27:01] 
    ## Drift correction results inspected, report:
    ## Drift_corrected: 100%

``` r
#corrected <- correct_drift(corrected)
# Flag low quality features
#corrected <- flag_quality(corrected, condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
```

Then, we can visualize the data after drift correction.

``` r
# Boxplot
corr_bp <- plot_sample_boxplots(corrected,
                                order_by = "Group",
                                fill_by = "Species")
# PCA
corr_pca <- plot_pca(corrected,
                     center = TRUE,
                     shape = "Species",
                     color = "Species")
# Package to plots visualization in a same windows
library(patchwork)
# Plot
corr_pca + corr_bp
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

The next step is feature clustering. This step helps us reduce the
number of features of the same molecule that were split due to a 70 eV
electron ionization (EI) environment.

``` r
clustered <- cluster_features(corrected, rt_window = 1/60, corr_thresh = 0.95,
                              d_thresh = 0.80)
```

    ## INFO [2025-10-13 08:27:03] Identified m/z column mass and retention time column RT
    ## INFO [2025-10-13 08:27:03] 
    ## Starting feature clustering at 2025-10-13 08:27:03.546195
    ## INFO [2025-10-13 08:27:03] Finding connections between features in RTX5MS_EI
    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500
    ## [1] 600
    ## INFO [2025-10-13 08:27:25] Found 5601 connections in RTX5MS_EI
    ## INFO [2025-10-13 08:27:25] Found 5601 connections
    ## 57 components found
    ## 
    ## 29 components found
    ## 
    ## 2 components found
    ## 
    ## 2 components found
    ## 
    ## 2 components found
    ## 
    ## INFO [2025-10-13 08:27:26] Found 72 clusters of 2 or more features, clustering finished at 2025-10-13 08:27:26.257494

``` r
compressed <- compress_clusters(clustered)
```

    ## INFO [2025-10-13 08:27:26] Clusters compressed, left with 179 features

We can inspect the data using the PCA plot after the clustering
algorithm execution.

``` r
# Boxplot
clust_bp <- plot_sample_boxplots(compressed,
                                 order_by = "Group",
                                 fill_by = "Species")
# PCA
clust_pca <- plot_pca(compressed,
                      center = TRUE,
                      shape = "Species",
                      color = "Species")
# Plot
clust_pca + clust_bp
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## Processing

We used probabilistic quotient normalization (PQN) with missing value
imputation and scaled using the autoscaling method for downstream
statistical analysis.

The code below imputes the data using a chromatographic noise value.

``` r
# Impute missing values using noise threshold
imputed <- impute_simple(compressed, value = 45, na_limit = 0)
```

We can inspect the data with the PCA plot after data imputation.

``` r
# Boxplot
imp_bp <- plot_sample_boxplots(imputed,
                               order_by = "Group",
                               fill_by = "Species")
# PCA
imp_pca <- plot_pca(imputed,
                    center = TRUE,
                    shape = "Species",
                    color = "Species")
# Plot
imp_pca + imp_bp
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

After data imputation, the data were normalized by PQN.

``` r
# Probabilistic quotient normalization
pqn_set <- pqn_normalization(imputed,
                             ref = c("qc", "all"),
                             method = c("median", "mean"),
                             all_features = FALSE)
```

    ## INFO [2025-10-13 08:27:28] Starting PQN normalization
    ## INFO [2025-10-13 08:27:28] Using median of qc samples as reference spectrum

We can inspect the data with the PCA plot after data normalization.

``` r
# Boxplot
pqn_bp <- plot_sample_boxplots(pqn_set,
                               order_by = "Group",
                               fill_by = "Species")
# PCA
pqn_pca <- plot_pca(pqn_set,
                    center = TRUE,
                    shape = "Species",
                    color = "Species")
# Plot
pqn_pca + pqn_bp
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Finally, the data is ready to be exported and proceed with the
statistical analysis.

``` r
#save(compressed, file = paste0(ppath, "Result/notame_Result/HS_GCMS/Notame_HS_GC-MS_out.RData"))
```

# Principal component analysis (PCA)

Extracting the data and metadata of “*notame*” MetaboSet for principal
component analysis (PCA) plot.

``` r
# Extract clean data
pqn_noflag <- drop_flagged(pqn_set)
# Extracting feature height table
peak_height <- exprs(pqn_noflag)
# Extracting phenotypic data
pheno_data <- pqn_noflag@phenoData@data
```

Transposing the feature table and preparing the PCA data.

``` r
# Transposing feature height table
transp_table  <- t(peak_height)
# Centering and Scaling features
ei_pca <- prcomp(transp_table, center = TRUE, scale. = TRUE)
```

## Score PCA

Plotting PCA results.

``` r
# Library to left_join use
library(dplyr)
# PCA scores
scores <- ei_pca$x %>%                   # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pheno_data)                  # Adding metadata
# PCA plot
pca_plot <- ggplot(scores,
       aes(PC1, PC2, shape = Species, color = Species)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (20.14 %)"),
         y=guide_axis(title = "PC2 (17.62 %)")) +
  labs(shape = 'Bees species', color= 'Bees species') +
  theme_classic() +
  theme(legend.text = element_text(face="italic")) +
  theme(legend.position = c(0.120, 0.230),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
pca_plot
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

## Loading PCA

Plotting a loading PCA result.

``` r
loadings <- ei_pca$rotation %>%           # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a column with
compound names (identified metabolites).

``` r
# Extracting feature identified
metab_data <- pqn_noflag[!is.na(pqn_noflag@featureData@data$Metabolite),]
# Extracting metabolite table
meta_table <- metab_data@featureData@data
# Creating a new small table of the annotated compounds
ei_compouds <- left_join(meta_table, loadings)
# Plotting results
load_pca <- ggplot(loadings, aes(PC1, PC2)) + 
  geom_point(alpha = 0.3, size = 2) +
  theme_classic() + 
  geom_point(data = ei_compouds,
             aes(shape = meta_table$IL,
                 color = meta_table$IL),
             size = 2.5) +
  labs(shape = 'Identification level',
       color = 'Identification level') +
  scale_color_manual(values = c("green",
                                "darkblue")) +
  scale_shape_manual(values = c(17, 19)) +
  ggrepel::geom_label_repel(data = ei_compouds,
                            aes(label = meta_table$Metabolite),
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 50,
                            min.segment.length = 0) +
  guides(x=guide_axis(title = "PC1 (20.14 %)"),
         y=guide_axis(title = "PC2 (17.62 %)")) +
  theme(legend.position = c(0.070, 0.090),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
  #ggsci::scale_color_aaas()
load_pca
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
# Save plot
#ggsave('Result/notame_Result/HS_GCMS/load_pca.pdf', width = 14, height = 7, device='pdf', dpi="print")
```

# PERMANOVA analysis

## *T. angustula*

PERMANOVA analysis of *T. angustula* by collection site.

``` r
# Drop QC
pqn_noqc <- drop_qcs(pqn_noflag)
# Removal of the M. fasciculata group from the dataset
permanova_ta <- pqn_noqc[, pqn_noqc$Species != "M. fasciculata"]
pData(permanova_ta) <- droplevels(pData(permanova_ta))
# Removal of the M. fuscopilosa group from the dataset
permanova_ta <- permanova_ta[, permanova_ta$Species != "M. fuscopilosa"]
pData(permanova_ta) <- droplevels(pData(permanova_ta))
# Change to factor varibles
permanova_ta@phenoData@data$Group <- as.factor(permanova_ta@phenoData@data$Group)
# PERMANOVA between close hives of bees T. angustula site 1 and 2 (Removal of site 3)
permanova_1vs2ta <- permanova_ta[, permanova_ta$Site != "3"]
pData(permanova_1vs2ta) <- droplevels(pData(permanova_1vs2ta))
permanova_1vs2ta <- perform_permanova(permanova_1vs2ta, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:33] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:34] PERMANOVA performed

``` r
permanova_1vs2ta
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           ER001 RYTA006
    ## C ER001       1       0
    ## C RYTA006     0       1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  384.4018 510.5982      1        4 3.011384   0.083        0.083

``` r
# PERMANOVA between close hives of bees T. angustula site 1 and 3 (Removal of site 2)
permanova_1vs3ta <- permanova_ta[, permanova_ta$Site != "2"]
pData(permanova_1vs3ta) <- droplevels(pData(permanova_1vs3ta))
permanova_1vs3ta <- perform_permanova(permanova_1vs3ta, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:34] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:34] PERMANOVA performed

``` r
permanova_1vs3ta
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           ATA001 RYTA006
    ## C ATA001       1       0
    ## C RYTA006      0       1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  543.1867 351.8133      1        4 6.175852     0.1          0.1

``` r
# PERMANOVA between close hives of bees T. angustula site 2 and 3 (Removal of site 1)
permanova_2vs3ta <- permanova_ta[, permanova_ta$Site != "1"]
pData(permanova_2vs3ta) <- droplevels(pData(permanova_2vs3ta))
permanova_2vs3ta <- perform_permanova(permanova_2vs3ta, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:34] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:35] PERMANOVA performed

``` r
permanova_2vs3ta
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##          ATA001 ER001
    ## C ATA001      1     0
    ## C ER001       0     1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  516.7884 378.2116      1        4 5.465602   0.082        0.082

## *M. fasciculata*

PERMANOVA analysis of *M. fasciculata* by collection site.

``` r
# Removal of the T. angustula group from the dataset
permanova_mfa <- pqn_noqc[, pqn_noqc$Species != "T. angustula"]
pData(permanova_mfa) <- droplevels(pData(permanova_mfa))
# Removal of the M. fuscopilosa group from the dataset
permanova_mfa <- permanova_mfa[, permanova_mfa$Species != "M. fuscopilosa"]
pData(permanova_mfa) <- droplevels(pData(permanova_mfa))
# Change to factor varibles
permanova_mfa@phenoData@data$Group <- as.factor(permanova_mfa@phenoData@data$Group)
# PERMANOVA between close hives of bees M. fasciculata site 5 and 6 (Removal of site 4)
permanova_5vs6mfa <- permanova_mfa[, permanova_mfa$Site != "4"]
pData(permanova_5vs6mfa) <- droplevels(pData(permanova_5vs6mfa))
permanova_5vs6mfa <- perform_permanova(permanova_5vs6mfa, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:35] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:36] PERMANOVA performed

``` r
permanova_5vs6mfa
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##          RGB004 RIG005
    ## C RGB004      1      0
    ## C RIG005      0      1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total   619.634  275.366      1        4 9.000877   0.094        0.094

``` r
# PERMANOVA between close hives of bees M. fasciculata site 5 and 4 (Removal of site 6)
permanova_5vs4mfa <- permanova_mfa[, permanova_mfa$Site != "6"]
pData(permanova_5vs4mfa) <- droplevels(pData(permanova_5vs4mfa))
permanova_5vs4mfa <- perform_permanova(permanova_5vs4mfa, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:36] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:36] PERMANOVA performed

``` r
permanova_5vs4mfa
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           RIG005 RYMG001
    ## C RIG005       1       0
    ## C RYMG001      0       1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  548.2971 346.7029      1        4 6.325843   0.078        0.078

``` r
# PERMANOVA between close hives of bees M. fasciculata site 6 and 4 (Removal of site 5)
permanova_6vs4mfa <- permanova_mfa[, permanova_mfa$Site != "5"]
pData(permanova_6vs4mfa) <- droplevels(pData(permanova_6vs4mfa))
permanova_6vs4mfa <- perform_permanova(permanova_6vs4mfa, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:36] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:37] PERMANOVA performed

``` r
permanova_6vs4mfa
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           RGB004 RYMG001
    ## C RGB004       1       0
    ## C RYMG001      0       1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  551.4406 343.5594      1        4 6.420323   0.093        0.093

## *M. fuscopilosa*

PERMANOVA analysis of *M. fuscopilosa* by collection site.

``` r
# Removal of the T. angustula group from the dataset
permanova_mfu <- pqn_noqc[, pqn_noqc$Species != "T. angustula"]
pData(permanova_mfu) <- droplevels(pData(permanova_mfu))
# Removal of the M. fasciculata group from the dataset
permanova_mfu <- permanova_mfu[, permanova_mfu$Species != "M. fasciculata"]
pData(permanova_mfu) <- droplevels(pData(permanova_mfu))
# Change to factor varibles
permanova_mfu@phenoData@data$Group <- as.factor(permanova_mfu@phenoData@data$Group)
# PERMANOVA between close hives of bees M. fuscopilosa site 8 and 9 (Removal of site 7)
permanova_8vs9mfu <- permanova_mfu[, permanova_mfu$Site != "7"]
pData(permanova_8vs9mfu) <- droplevels(pData(permanova_8vs9mfu))
permanova_8vs9mfu <- perform_permanova(permanova_8vs9mfu, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:37] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:38] PERMANOVA performed

``` r
permanova_8vs9mfu
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##          RGY001 RIG003
    ## C RGY001      1      0
    ## C RIG003      0      1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  554.0464 340.9536      1        4 6.499961   0.081        0.081

``` r
# PERMANOVA between close hives of bees M. fuscopilosa site 8 and 7 (Removal of site 9)
permanova_8vs7mfu <- permanova_mfu[, permanova_mfu$Site != "9"]
pData(permanova_8vs7mfu) <- droplevels(pData(permanova_8vs7mfu))
permanova_8vs7mfu <- perform_permanova(permanova_8vs7mfu, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:38] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:38] PERMANOVA performed

``` r
permanova_8vs7mfu
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           FIED001 RGY001
    ## C FIED001       1      0
    ## C RGY001        0      1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  531.4305 363.5695      1        4 5.846809   0.046        0.046

``` r
# PERMANOVA between close hives of bees M. fuscopilosa site 9 and 7 (Removal of site 8)
permanova_9vs7mfu <- permanova_mfu[, permanova_mfu$Site != "8"]
pData(permanova_9vs7mfu) <- droplevels(pData(permanova_9vs7mfu))
permanova_9vs7mfu <- perform_permanova(permanova_9vs7mfu, group = "Group", nperm = 999)
```

    ## INFO [2025-10-13 08:27:38] Starting PERMANOVA tests
    ## INFO [2025-10-13 08:27:39] PERMANOVA performed

``` r
permanova_9vs7mfu
```

    ##  ###### PERMANOVA Analysis #######
    ## 
    ## Call
    ## PERMANOVA::PERMANOVA(Distance = initialized, group = pData(object)[, 
    ##     group], nperm = 999)
    ## ________________________________________________
    ## 
    ## Contrast Matrix
    ##           FIED001 RIG003
    ## C FIED001       1      0
    ## C RIG003        0      1
    ## ________________________________________________
    ## 
    ## MANOVA
    ##       Explained Residual df Num df Denom    F-exp p-value p-value adj.
    ## Total  607.7088 287.2912      1        4 8.461223     0.1          0.1

# Heatmap with HCA

Only the annotated features were used in the heatmap and the
hierarchical cluster analysis (HCA). The ComplexHeatmap package will be
used for the heatmap and the hierarchical cluster analysis (HCA).

ComplexHeatmap package and dependency installation.

``` r
# ComplexHeatmap package installation
#if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

# ColorRamp2 package installation
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("jokergoo/colorRamp2")
library(colorRamp2)

# Cowplot package installation
#install.packages("cowplot")
library(cowplot)

# mdatools package installation
#install_github('svkucheryavski/mdatools')
library(mdatools)

# ClassyFire package installation
#remotes::install_github('aberHRML/classyfireR')
library(classyfireR)
```

The metabolites were classified using the
[ClassyFireR](https://doi.org/10.1186/s13321-016-0174-y) to add these
metabolite classifications to the heatmap plot.

``` r
# InChI key of the metabolites you want to classify
InChI_Keys <- c('2-Heptanone' = "CATSNJVOTSVZJV-UHFFFAOYSA-N")
# Get classification
Classification_List <- purrr::map(InChI_Keys, get_classification)
Classification_List
```

Extracting and loading of identified metabolites’ abundance, and data
scaling by the autoscaling method.

``` r
# Drop QC
hm_no_qc <- drop_qcs(pqn_noflag)
# Scaling by autoscaling method
hm_scl <- scale(t(exprs(hm_no_qc)), center = TRUE, scale = TRUE)
hm_scl <- t(hm_scl)
# Adding autoscaling data to notame MetaboSet
hm_scl_set <- hm_no_qc
exprs(hm_scl_set) <- hm_scl
# Extracting identified metabolite data
raw_hm <- hm_scl_set[!is.na(hm_scl_set@featureData@data$Metabolite),]
# Extracting feature height table
hm_height <- exprs(raw_hm)
# Extracting sample information
hm_pdata <- raw_hm@phenoData@data
# Extracting feature information
hm_fdata <- raw_hm@featureData@data
```

Row and top heatmap annotation.

``` r
set.seed(1540)
# Adding row and column names
hm_scl <- hm_height
rownames(hm_scl) <- hm_fdata$Metabolite
colnames(hm_scl) <- hm_pdata$Species
# Metabolite class color
cols_metclass <- c("Benzenoids" = "#800000FF",
                   "Hydrocarbons" = "#FFA319FF",
                   "Lipids and lipid-like molecules" = "#8A9045FF",
                   "Organic oxygen compounds" = "#8DD3C7",
                   "Organohalogen compounds" = "#BEBADA",
                   "Organoheterocyclic compounds" = "#FFFFB3")
# Add row anotation to HeatMap
hm_row_ann <- rowAnnotation(`Superclass` = hm_fdata$classyfireR_Superclass,
                            col = list(`Superclass` = cols_metclass),
                            show_annotation_name = T,
                            show_legend = F)
# Species color
cols_species <- c("T. angustula" = "#E76BF3",
                  "M. fasciculata" = "#F8766D",
                  "M. fuscopilosa" = "#7CAE00")
# Add top anotation to Heatmap
top_info_ann <- HeatmapAnnotation(`Species` = hm_pdata$Species,
                                  col = list(`Species` = cols_species),
                                  show_annotation_name = T,
                                  show_legend = F,
                                  border = T)
# Color scale
mycol <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
# Heatmap matrix plotting
hm_plot <- Heatmap(hm_scl,
                   col = mycol,
                   border_gp = grid::gpar(col = "black", lty = 0.02),
                   rect_gp = grid::gpar(col = "black", lwd = 0.75),
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "complete",
                   top_annotation = top_info_ann,
                   column_names_gp = gpar(fontface = "italic"),
                   row_names_max_width = unit(10, "cm"),
                   right_annotation = hm_row_ann,
                   show_heatmap_legend = F,
                   row_km = 3, column_km = 2,
                   row_title = c("a", "b", "c"))
hm_plot
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

Adding legends to the heatmap.

``` r
# Color scale legend
lgd1 <- Legend(col_fun = mycol,
               title = "Autoscaled abundance",
               direction = "horizontal" )
# Bees species legend
lgd2 <- Legend(labels = gt_render(c("*T. angustula*",
                                    "*M. fasciculata*",
                                    "*M. fuscopilosa*")),
               legend_gp = gpar(fill = cols_species),
               title = "Bees species", ncol = 1)
# Metabolite class Legend
lgd3 <- Legend(labels = c(unique(hm_fdata$classyfireR_Superclass)) ,
               legend_gp = gpar(fill = cols_metclass), 
               title = "Metabolite superclass", ncol = 2)
```

Plotting the heatmap.

``` r
set.seed(1540)
# Converting to ggplot
gg_heatmap <- grid.grabExpr(draw(hm_plot))
gg_heatmap <- ggpubr::as_ggplot(gg_heatmap)
# Legends
all_legends <- packLegend(lgd1, lgd2, lgd3, direction = "horizontal")
gg_legend <- grid.grabExpr(draw(all_legends))
gg_legend_fn <- ggpubr::as_ggplot(gg_legend)
# Heatmap plot
gcms_hm <- plot_grid(gg_legend_fn,
                     gg_heatmap, ncol = 1,
                     rel_heights = c(0.055, 0.880))
gcms_hm
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

# Antimicrobial activity correlation

Correlation between feature abundance and antibacterial result (zone
inhibition diameter) using the PCA and Pearson’s correlation approaches.

## Using the PCA

Data and metadata preparation, merging the antimicrobial data with
feature abundance.

``` r
# Antimicrobial activity data extraction
bacteria_data <- drop_qcs(modes$na_na)
bacteria_data <- t(exprs(bacteria_data))
# Drop QCs of metabolomics data
no_qc <- drop_qcs(pqn_noflag)
# Extracting feature height table
corr_peak <- t(exprs(no_qc))
# Merge the antimicrobial activity result with the feature list table
bacteria_peak <- cbind(bacteria_data, corr_peak)
# Extracting phenotypic data of metabolomics data
corr_pdata <- no_qc@phenoData@data
```

### Score PCA

PCA scores calculation of the antimicrobial and feature abundance merge
data.

``` r
# Centering and Scaling features
corr_pca <- prcomp(bacteria_peak, center = TRUE, scale. = TRUE)
```

Plotting PCA score.

``` r
# PCA scores
corr_scores <- corr_pca$x %>%            # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(corr_pdata)                  # Adding metadata
# PCA plot
corr_pca_plot <- ggplot(corr_scores,
                        aes(PC1, PC2, shape = Species, color = Species)) +
  scale_color_manual(values=c("#F8766D",
                              "#7CAE00",
                              "#E76BF3")) +
  scale_shape_manual(values=c(16, 17, 3)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (22.46 %)"),
         y=guide_axis(title = "PC2 (18.82 %)")) +
  labs(shape = 'Bees species', color= 'Bees species') +
  theme_classic() +
  theme(legend.text = element_text(face="italic")) +
  theme(legend.position = c(0.120, 0.200),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
corr_pca_plot
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

### Loading PCA

Plotting loading results.

``` r
corr_loadings <- corr_pca$rotation %>%    # Extract loadings
  data.frame(Feature_ID = rownames(.))    # New column with feat name
```

Creating an artificial table with feature names and a column with
compound names (identified metabolites).

``` r
# Extracting the antimicrobial metadata
bacteria_table <- drop_qcs(modes$na_na)@featureData@data
bacteria_table$Flag  <- as.character(bacteria_table$Flag)
# Merge the antimicrobial metadata with the metabolite table
corr_table <- merge(meta_table, bacteria_table, all = TRUE)
# Creating a new small table of the annotated compounds and antimicrobial metadata
corr_compouds <- left_join(corr_table, corr_loadings)
# Plotting results
corr_load_plot <- ggplot(corr_loadings, aes(PC1, PC2)) + 
  geom_point(alpha = 0.3, size = 2) +
  theme_classic() +
  geom_point(data = corr_compouds, aes(shape = corr_table$IL,
                                       color = corr_table$IL), size = 2.5) +
  labs(shape = 'Legend of features',
       color = 'Legend of features') +
  scale_color_manual(values = c("red",
                                "green",
                                "darkblue")) +
  scale_shape_manual(values = c(15, 17, 19)) +
  ggrepel::geom_label_repel(data = corr_compouds,
                            aes(label = corr_table$Metabolite,
                                fontface = ifelse(corr_table$IL == "Antimicrobial activity",
                                                  'italic', 'plain')),
                            box.padding = 0.37,
                            label.padding = 0.22,
                            label.r = 0.30,
                            cex = 2.5,
                            max.overlaps = 50,
                            min.segment.length = 0) +
  guides(x=guide_axis(title = "PC1 (22.46 %)"),
         y=guide_axis(title = "PC2 (18.82 %)")) +
  theme(legend.position = c(0.070, 0.107),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
  #ggsci::scale_color_aaas()
corr_load_plot
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

## Using the Pearson correlation

The Pearson correlation between antimicrobial activity and feature
abundance was performed.

### Normality of the antimicrobial data

First, the antimicrobial data normality was evaluated using the
Shapiro-Wilk test.

``` r
# Loading the "readxl" package
library(readxl)
# Loading the antimicrobial activity data
all_bacteria <- read_excel("../Data/Antimicrobial_of_honey.xlsx", sheet = 1)
# Deleting the group column to perform Shapiro-Wilk test in batch
no_group <- all_bacteria[,-1]
# Perform the Shapiro-Wilk test
shapiro <- do.call(rbind, lapply(no_group,
                                 function(x) shapiro.test(x)[c("statistic", "p.value")]))
shapiro
```

    ##                 statistic p.value     
    ## E_coli_100      0.8894981 0.007754756 
    ## E_coli_75       0.910591  0.0235632   
    ## K_pneumonia_50  0.9506033 0.2217358   
    ## K_pneumonia_75  0.9520331 0.2400907   
    ## K_pneumonia_100 0.8965943 0.01118818  
    ## P_mirabilis_50  0.9317354 0.07623142  
    ## P_mirabilis_75  0.9404336 0.1248864   
    ## P_mirabilis_100 0.9253929 0.05333025  
    ## S_aureus_25     0.7805308 6.412762e-05
    ## S_aureus_50     0.8894501 0.007735745 
    ## S_aureus_75     0.8934946 0.009524042 
    ## S_aureus_100    0.896372  0.01105921

The result showed that only the microorganisms *K. pneumonia* and *P.
mirabilis* have a normal distribution (p-value \> 0.05). However, when
we inspected the data, we found the values in the ATA001 and RIG003
groups in the *E. coli* could be outliers; also, a possible outlier in
the RYMG001 group in the *S. aureus* was found. When we delete all these
possible outliers, we find that all the microorganisms (not including
*S. aureus* at 25 % v/v) have a normal distribution, as seen below.

``` r
# Delete outlier data in E. coli
no_outliers <- all_bacteria
no_outliers[no_outliers$Group == "ATA001",]$E_coli_75 <- NA
no_outliers[no_outliers$Group == "ATA001",]$E_coli_100 <- NA
no_outliers[no_outliers$Group == "RIG003",]$E_coli_75 <- NA
no_outliers[no_outliers$Group == "RIG003",]$E_coli_100 <- NA
# Delete outlier data in S. aureus
no_outliers[no_outliers$Group == "RYMG001",]$S_aureus_25 <- NA
no_outliers[no_outliers$Group == "RYMG001",]$S_aureus_50 <- NA
no_outliers[no_outliers$Group == "RYMG001",]$S_aureus_75 <- NA
no_outliers[no_outliers$Group == "RYMG001",]$S_aureus_100 <- NA
# Deleting the group column to perform Shapiro-Wilk test in batch
no_outliers <- no_outliers[,-1]
# Perform the Shapiro-Wilk test
shapiro_outl <- do.call(rbind, lapply(no_outliers,
                                 function(x) shapiro.test(x)[c("statistic", "p.value")]))
shapiro_outl
```

    ##                 statistic p.value    
    ## E_coli_100      0.9319303 0.1504756  
    ## E_coli_75       0.9202502 0.08778051 
    ## K_pneumonia_50  0.9506033 0.2217358  
    ## K_pneumonia_75  0.9520331 0.2400907  
    ## K_pneumonia_100 0.8965943 0.01118818 
    ## P_mirabilis_50  0.9317354 0.07623142 
    ## P_mirabilis_75  0.9404336 0.1248864  
    ## P_mirabilis_100 0.9253929 0.05333025 
    ## S_aureus_25     0.8831963 0.009642998
    ## S_aureus_50     0.9319041 0.1075295  
    ## S_aureus_75     0.9669953 0.5936641  
    ## S_aureus_100    0.9608898 0.4565626

The antimicrobial activity result of honey at 75 % v/v was used in the
further analysis because this concentration showed activity in most of
the tested microorganisms. (*E. coli*, *K. pneumonia*, *P. mirabilis*,
and *S. aureus*).

### Metabolomics data transformation

Before the Pearson correlation of feature abundance with the
antimicrobial activity (zone inhibition diameter), the previously
normalized data by PQN was transformed using the generalised logarithm
(glog) method.

``` r
# "SummarizedExperiment" package installation
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
# Convert feature height table to SummarizedExperiment class
pmp_data <- SummarizedExperiment(assays = exprs(pqn_noflag),
                                 colData = pqn_noflag@phenoData@data)
# Package for generalized logarithmic transform
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("pmp")
library(pmp)
# Generalised logarithmic transform
glog_exprs <- glog_transformation(df = pmp_data@assays@data@listData[[1]],
                                  classes = pmp_data$QC,
                                  qc_label = "QC")
# Adding glog transformation to notame MetaboSet
glog_set <- pqn_noflag
exprs(glog_set) <- glog_exprs
```

Preparing the metabolomics data for the antimicrobial correlation.

``` r
# Drop QCs
glog_no_qc <- drop_qcs(glog_set)
# Preparing the data for E. coli correlation (drop outlier data)
glog_E.coli <- glog_no_qc[, glog_no_qc$Group != "ATA001"]
pData(glog_E.coli) <- droplevels(pData(glog_E.coli))
glog_E.coli <- glog_E.coli[, glog_E.coli$Group != "RIG003"]
pData(glog_E.coli) <- droplevels(pData(glog_E.coli))
# Preparing the data for S. aureus (dropn outlier data)
glog_S.aureus <- glog_no_qc[, glog_no_qc$Group != "RYMG001"]
pData(glog_S.aureus) <- droplevels(pData(glog_S.aureus))
# Convert antimicrobial result to numeric
glog_E.coli$E_coli_75 <- as.numeric(glog_E.coli$E_coli_75)
glog_no_qc$K_pneumonia_75 <- as.numeric(glog_no_qc$K_pneumonia_75)
glog_no_qc$P_mirabilis_75 <- as.numeric(glog_no_qc$P_mirabilis_75)
glog_S.aureus$S_aureus_75 <- as.numeric(glog_S.aureus$S_aureus_75)
```

Perform the Pearson correlation analysis.

``` r
# Feature correlation with the E. coli antimicrobial activity
E.coli_corr <- perform_correlation_tests(glog_E.coli,
                                         x = featureNames(glog_no_qc),
                                         y = c("E_coli_75"),
                                         method = "pearson")
```

    ## INFO [2025-10-13 08:27:53] Starting correlation tests.
    ## INFO [2025-10-13 08:27:53] Correlation tests performed.

``` r
# Feature correlation with the K. pneumonia antimicrobial activity
K.pneumonia_corr <- perform_correlation_tests(glog_no_qc,
                                              x = featureNames(glog_set),
                                              y = c("K_pneumonia_75"),
                                              method = "pearson")
```

    ## INFO [2025-10-13 08:27:53] Starting correlation tests.
    ## INFO [2025-10-13 08:27:54] Correlation tests performed.

``` r
# Feature correlation with the P. mirabilis antimicrobial activity
P.mirabilis_corr <- perform_correlation_tests(glog_no_qc,
                                              x = featureNames(glog_set),
                                              y = c("P_mirabilis_75"),
                                              method = "pearson")
```

    ## INFO [2025-10-13 08:27:54] Starting correlation tests.
    ## INFO [2025-10-13 08:27:54] Correlation tests performed.

``` r
# Feature correlation with the S. aureus antimicrobial activity
S.aureus_corr <- perform_correlation_tests(glog_S.aureus,
                                           x = featureNames(glog_no_qc),
                                           y = c("S_aureus_75"),
                                           method = "pearson")
```

    ## INFO [2025-10-13 08:27:54] Starting correlation tests.
    ## INFO [2025-10-13 08:27:54] Correlation tests performed.

### Correlation of identified metabolites

Create a matrix for the heatplot of Pearson’s correlation.

``` r
# Correlation table for the heatmap
hm_corr_table <- data.frame(Feature_ID = E.coli_corr$X,
                            E_coli = E.coli_corr$Correlation_coefficient,
                            K_pneumonia = K.pneumonia_corr$Correlation_coefficient,
                            P_mirabilis = P.mirabilis_corr$Correlation_coefficient,
                            S_aureus = S.aureus_corr$Correlation_coefficient)
# Feature data table
feat_corr_table <- data.frame(Feature_ID = hm_fdata$Feature_ID,
                             Metabolite = hm_fdata$Metabolite)
# Adding the metabolie name to the correlation table
hm_corr_table <- left_join(hm_corr_table, feat_corr_table)
# Extracting data of the identified metabolites
hm_corr_table <- hm_corr_table[!is.na(hm_corr_table$Metabolite),]
# Adding row name
rownames(hm_corr_table) <- hm_corr_table$Metabolite
# Delete extra information
hm_corr_table <- subset(hm_corr_table, select = -c(Feature_ID, Metabolite))
# Converting DataFrame to data matrix
hm_corr_table <- data.matrix(hm_corr_table, rownames.force = NA)
# Adding column name
colnames(hm_corr_table) <- c("E. coli", "K. pneumonia", "P. mirabilis", "S. aureus")
```

Create a matrix showing the statistical significance of the correlation
result.

``` r
# Correlation table for the heatmap
hm_sign_table <- data.frame(Feature_ID = E.coli_corr$X,
                            E_coli = E.coli_corr$Correlation_P_FDR,
                            K_pneumonia = K.pneumonia_corr$Correlation_P_FDR,
                            P_mirabilis = P.mirabilis_corr$Correlation_P_FDR,
                            S_aureus = S.aureus_corr$Correlation_P_FDR)
# Adding the metabolie name to the correlation table
hm_sign_table <- left_join(hm_sign_table, feat_corr_table)
# Extracting data of the identified metabolites
hm_sign_table <- hm_sign_table[!is.na(hm_sign_table$Metabolite),]
# Adding row name
rownames(hm_sign_table) <- hm_sign_table$Metabolite
# Delete extra information
hm_sign_table <- subset(hm_sign_table, select = -c(Feature_ID, Metabolite))
# Converting DataFrame to data matrix
hm_sign_table <- data.matrix(hm_sign_table, rownames.force = NA)
# Adding column name
colnames(hm_sign_table) <- c("E. coli", "K. pneumonia", "P. mirabilis", "S. aureus")
```

Plot the correlation heatplot.

``` r
# Add top anotation to Heatmap
ann_corr <- HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = 0, col = "white"),
                   labels = "Pearson correlation",
                   labels_gp = gpar(col = "black", fontsize = 12,
                                    fontface = "bold")))
# Color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# Heatmap plot
hm_corr <- Heatmap(hm_corr_table, col = col_fun, cluster_rows = FALSE,
                   cluster_columns = FALSE, column_names_gp = gpar(fontface = "italic"),
                   right_annotation = hm_row_ann,
                   border_gp = grid::gpar(col = "black", lty = 0.02),
                   rect_gp = grid::gpar(col = "black", lwd = 0.75),
                   top_annotation = ann_corr,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(hm_sign_table[i, j] < 0.05) {
                       grid.text(sprintf("%.1f", hm_corr_table[i, j]), x, y,
                                 gp = gpar(fontsize = 12, fontface = "bold"))
                       } else if(hm_sign_table[i, j] > 0.05) {
                         grid.text(sprintf("%.1f", hm_corr_table[i, j]), x, y,
                                   gp = gpar(fontsize = 12))}},
                   show_heatmap_legend = F)
hm_corr
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

### Correlation of unknown metabolites

Filtering the metabolites with positive and high Pearson’s correlation
coefficient (r \>0.5).

``` r
# Library to change character values
library('stringr')
# Correlation table for the heatmap with the best features
best_feat <- data.frame(Feature_ID = E.coli_corr$X,
                        E_coli = E.coli_corr$Correlation_coefficient,
                        E_coli_p = E.coli_corr$Correlation_P_FDR,
                        K_pneumonia = K.pneumonia_corr$Correlation_coefficient,
                        K_pneumonia_p = K.pneumonia_corr$Correlation_P_FDR,
                        P_mirabilis = P.mirabilis_corr$Correlation_coefficient,
                        P_mirabilis_p = P.mirabilis_corr$Correlation_P_FDR,
                        S_aureus = S.aureus_corr$Correlation_coefficient,
                        S_aureus_p = S.aureus_corr$Correlation_P_FDR)
# Filtering the feature with the correlation >= 0.5
best_feat <- best_feat[best_feat$E_coli >= 0.5 |
                         best_feat$K_pneumonia >= 0.5 |
                         best_feat$P_mirabilis >= 0.5 |
                         best_feat$S_aureus >= 0.5, ]
# Adding the metabolite name to the correlation table
best_feat <- left_join(best_feat, feat_corr_table)
# Extracting data of no identified metabolites
best_feat <- best_feat[!complete.cases(best_feat$Metabolite),]
# Replace the column name and ionization type "RTX5MS_EI" with unknown
best_feat$Feature_ID <- sub('RTX5MS_EI_','Unknown-', best_feat$Feature_ID)
best_feat$Feature_ID <- str_replace_all(best_feat$Feature_ID, '_', '.')
```

Create a matrix for the heatplot of Pearson’s correlation.

``` r
# Correlation table for the heatmap
hm_corr_table1 <- data.frame(Feature_ID = best_feat$Feature_ID,
                             E_coli = best_feat$E_coli,
                             K_pneumonia = best_feat$K_pneumonia,
                             P_mirabilis = best_feat$P_mirabilis,
                             S_aureus = best_feat$S_aureus)
# Adding row name
rownames(hm_corr_table1) <- hm_corr_table1$Feature_ID
# Delete extra information
hm_corr_table1 <- subset(hm_corr_table1, select = -c(Feature_ID))
# Converting DataFrame to data matrix
hm_corr_table1 <- data.matrix(hm_corr_table1, rownames.force = NA)
# Adding column name
colnames(hm_corr_table1) <- c("E. coli", "K. pneumonia", "P. mirabilis", "S. aureus")
```

Create a matrix showing the statistical significance of the correlation
result.

``` r
# Correlation table for the heatmap
hm_sign_table1 <- data.frame(Feature_ID = best_feat$Feature_ID,
                             E_coli = best_feat$E_coli_p,
                             K_pneumonia = best_feat$K_pneumonia_p,
                             P_mirabilis = best_feat$P_mirabilis_p,
                             S_aureus = best_feat$S_aureus_p)
# Adding row name
rownames(hm_sign_table1) <- hm_sign_table1$Feature_ID
# Delete extra information
hm_sign_table1 <- subset(hm_sign_table1, select = -c(Feature_ID))
# Converting DataFrame to data matrix
hm_sign_table1 <- data.matrix(hm_sign_table1, rownames.force = NA)
# Adding column name
colnames(hm_sign_table1) <- c("E. coli", "K. pneumonia", "P. mirabilis", "S. aureus")
```

Plot the correlation heatplot.

``` r
# Color scale legend
lgd1a <- Legend(col_fun = col_fun,
                title = "Correlation coefficient",
                direction = "horizontal" )
# Packed legends
pd = packLegend(lgd1a, direction = "horizontal")
# Heatmap plot
#png(file="Result/notame_Result/HS_GCMS/FigureA1.png", width = 5, height = 6, units = "in", res= 300)
hm_corr1 <- Heatmap(hm_corr_table1, col = col_fun, cluster_rows = FALSE,
                   cluster_columns = FALSE, column_names_gp = gpar(fontface = "italic"),
                   border_gp = grid::gpar(col = "black", lty = 0.02),
                   rect_gp = grid::gpar(col = "black", lwd = 0.75),
                   top_annotation = ann_corr,
                   cell_fun = function(j, i, x, y, w, h, fill) {
                     if(hm_sign_table1[i, j] < 0.05) {
                       grid.text(sprintf("%.1f", hm_corr_table1[i, j]), x, y,
                                 gp = gpar(fontsize = 10, fontface = "bold"))
                       } else if(hm_sign_table1[i, j] > 0.05) {
                         grid.text(sprintf("%.1f", hm_corr_table1[i, j]), x, y,
                                   gp = gpar(fontsize = 10))}},
                   show_heatmap_legend = F)
hm_corr1
draw(pd, x = unit(7.5, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

``` r
#dev.off()
```

# UpSet plot

The UpSet plot was implemented to inspect the unique or shared features
of the honeys from each bee species and collection site (beehive).

``` r
# Mean of technical replicate
upset_aver <- summary_statistics(compressed, grouping_cols = "Group")
# Dataframe with average values by the collection sites (beehive)
upset_mat <- data.frame(Feature_ID = upset_aver$Feature_ID,
                        ata001 = upset_aver$ATA001_mean,
                        er001 = upset_aver$ER001_mean,
                        fied001 = upset_aver$FIED001_mean,
                        rgb004 = upset_aver$RGB004_mean,
                        rgy001 = upset_aver$RGY001_mean,
                        rig003 = upset_aver$RIG003_mean,
                        rig005 = upset_aver$RIG005_mean,
                        rymg001 = upset_aver$RYMG001_mean,
                        ryta006 = upset_aver$RYTA006_mean)
# Adding rownames
row.names(upset_mat) <- upset_mat$Feature_ID
upset_mat$Feature_ID <- NULL
# Adding colnames
colnames(upset_mat) <- c("T. angustula ATA001",
                         "T. angustula ER001",
                         "M. fuscopilosa FIED001",
                         "M. fasciculata RGB004",
                         "M. fuscopilosa RGY001",
                         "M. fuscopilosa RIG003",
                         "M. fasciculata RIG005",
                         "M. fasciculata RYMG001",
                         "T. angustula RYTA006")
# Change the values by presence using one (1) or absence using 0
upset_mat[upset_mat > 0] <- 1       # presence 
upset_mat[is.na(upset_mat)] <- 0    # absence
# Make the combination matrix
comb_mat = make_comb_mat(upset_mat)
#png(file="Result/notame_Result/HS_GCMS/FigureA2_1.png", width = 7, height = 3.5, units = "in", res= 300)
upset_plot <- UpSet(comb_mat,
                    row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                    top_annotation =
                      upset_top_annotation(comb_mat,
                                           gp = gpar(fill = "#F8766D",
                                                     col = "#F8766D"),
                                           add_numbers = TRUE,
                                           annotation_name_rot = 90,),
                    right_annotation =
                      upset_right_annotation(comb_mat,
                                             gp = gpar(fill = "#1C6AA8",
                                                     col = "#1C6AA8"),
                                             add_numbers = TRUE),
                    set_order = c("T. angustula ATA001",
                                  "T. angustula ER001",
                                  "T. angustula RYTA006",
                                  "M. fuscopilosa FIED001",
                                  "M. fuscopilosa RGY001",
                                  "M. fuscopilosa RIG003",
                                  "M. fasciculata RGB004",
                                  "M. fasciculata RIG005",
                                  "M. fasciculata RYMG001"))
#pdf(file="Result/notame_Result/HS_GCMS/FigureA2.pdf", width = 7, height = 3.5)
upset_plot
```

![](HS_GCMS_Multivariate_Statistics_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

``` r
#dev.off()
```

Finish a record.

``` r
finish_log()
```

    ## INFO [2025-10-13 08:27:57] Finished analysis. Mon Oct 13 08:27:57 2025
    ## Session info:
    ## 
    ## INFO [2025-10-13 08:27:57] R version 4.4.1 (2024-06-14 ucrt)
    ## INFO [2025-10-13 08:27:57] Platform: x86_64-w64-mingw32/x64
    ## INFO [2025-10-13 08:27:57] Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] Matrix products: default
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] locale:
    ## INFO [2025-10-13 08:27:57] [1] LC_COLLATE=English_United States.1252 
    ## INFO [2025-10-13 08:27:57] [2] LC_CTYPE=English_United States.1252   
    ## INFO [2025-10-13 08:27:57] [3] LC_MONETARY=English_United States.1252
    ## INFO [2025-10-13 08:27:57] [4] LC_NUMERIC=C                          
    ## INFO [2025-10-13 08:27:57] [5] LC_TIME=English_United States.1252    
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] time zone: America/Guayaquil
    ## INFO [2025-10-13 08:27:57] tzcode source: internal
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] attached base packages:
    ## INFO [2025-10-13 08:27:57] [1] stats4    grid      stats     graphics  grDevices utils     datasets 
    ## INFO [2025-10-13 08:27:57] [8] methods   base     
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] other attached packages:
    ## INFO [2025-10-13 08:27:57]  [1] stringr_1.5.1               pmp_1.14.1                 
    ## INFO [2025-10-13 08:27:57]  [3] SummarizedExperiment_1.34.0 GenomicRanges_1.56.0       
    ## INFO [2025-10-13 08:27:57]  [5] GenomeInfoDb_1.40.1         IRanges_2.38.0             
    ## INFO [2025-10-13 08:27:57]  [7] S4Vectors_0.42.0            MatrixGenerics_1.16.0      
    ## INFO [2025-10-13 08:27:57]  [9] matrixStats_1.3.0           readxl_1.4.3               
    ## INFO [2025-10-13 08:27:57] [11] classyfireR_0.3.8           mdatools_0.14.2            
    ## INFO [2025-10-13 08:27:57] [13] cowplot_1.1.3               colorRamp2_0.1.0           
    ## INFO [2025-10-13 08:27:57] [15] ComplexHeatmap_2.20.0       dplyr_1.1.4                
    ## INFO [2025-10-13 08:27:57] [17] patchwork_1.3.0             notame_0.3.1               
    ## INFO [2025-10-13 08:27:57] [19] magrittr_2.0.3              ggplot2_3.5.1              
    ## INFO [2025-10-13 08:27:57] [21] futile.logger_1.4.3         Biobase_2.64.0             
    ## INFO [2025-10-13 08:27:57] [23] BiocGenerics_0.54.0         generics_0.1.3             
    ## INFO [2025-10-13 08:27:57] 
    ## INFO [2025-10-13 08:27:57] loaded via a namespace (and not attached):
    ## INFO [2025-10-13 08:27:57]   [1] RColorBrewer_1.1-3      sys_3.4.2               rstudioapi_0.16.0      
    ## INFO [2025-10-13 08:27:57]   [4] jsonlite_1.8.8          shape_1.4.6.1           farver_2.1.2           
    ## INFO [2025-10-13 08:27:57]   [7] rmarkdown_2.27          GlobalOptions_0.1.2     fs_1.6.4               
    ## INFO [2025-10-13 08:27:57]  [10] zlibbioc_1.50.0         vctrs_0.6.5             memoise_2.0.1          
    ## INFO [2025-10-13 08:27:57]  [13] askpass_1.2.0           rstatix_0.7.2           itertools_0.1-3        
    ## INFO [2025-10-13 08:27:57]  [16] htmltools_0.5.8.1       S4Arrays_1.4.1          usethis_2.2.3          
    ## INFO [2025-10-13 08:27:57]  [19] missForest_1.5          lambda.r_1.2.4          curl_5.2.1             
    ## INFO [2025-10-13 08:27:57]  [22] broom_1.0.6             cellranger_1.1.0        SparseArray_1.4.8      
    ## INFO [2025-10-13 08:27:57]  [25] plyr_1.8.9              impute_1.80.0           futile.options_1.0.1   
    ## INFO [2025-10-13 08:27:57]  [28] cachem_1.1.0            commonmark_1.9.1        igraph_2.0.3           
    ## INFO [2025-10-13 08:27:57]  [31] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
    ## INFO [2025-10-13 08:27:57]  [34] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
    ## INFO [2025-10-13 08:27:57]  [37] GenomeInfoDbData_1.2.12 clue_0.3-65             digest_0.6.36          
    ## INFO [2025-10-13 08:27:57]  [40] pcaMethods_1.96.0       colorspace_2.1-0        RSQLite_2.3.7          
    ## INFO [2025-10-13 08:27:57]  [43] ggpubr_0.6.0            labeling_0.4.3          randomForest_4.7-1.1   
    ## INFO [2025-10-13 08:27:57]  [46] fansi_1.0.6             httr_1.4.7              abind_1.4-5            
    ## INFO [2025-10-13 08:27:57]  [49] compiler_4.4.1          rngtools_1.5.2          bit64_4.0.5            
    ## INFO [2025-10-13 08:27:57]  [52] withr_3.0.0             doParallel_1.0.17       backports_1.5.0        
    ## INFO [2025-10-13 08:27:57]  [55] carData_3.0-5           DBI_1.2.3               highr_0.11             
    ## INFO [2025-10-13 08:27:57]  [58] ggsignif_0.6.4          MASS_7.3-60.2           openssl_2.2.0          
    ## INFO [2025-10-13 08:27:57]  [61] DelayedArray_0.30.1     rjson_0.2.21            tools_4.4.1            
    ## INFO [2025-10-13 08:27:57]  [64] zip_2.3.1               glue_1.8.0              gridtext_0.1.5         
    ## INFO [2025-10-13 08:27:57]  [67] reshape2_1.4.4          cluster_2.1.6           gtable_0.3.5           
    ## INFO [2025-10-13 08:27:57]  [70] tidyr_1.3.1             xml2_1.3.6              car_3.1-2              
    ## INFO [2025-10-13 08:27:57]  [73] utf8_1.2.4              XVector_0.44.0          ggrepel_0.9.5          
    ## INFO [2025-10-13 08:27:57]  [76] foreach_1.5.2           pillar_1.9.0            markdown_1.13          
    ## INFO [2025-10-13 08:27:57]  [79] circlize_0.4.16         lattice_0.22-6          bit_4.0.5              
    ## INFO [2025-10-13 08:27:57]  [82] deldir_2.0-4            tidyselect_1.2.1        knitr_1.47             
    ## INFO [2025-10-13 08:27:57]  [85] xfun_0.45               credentials_2.0.1       stringi_1.8.4          
    ## INFO [2025-10-13 08:27:57]  [88] UCSC.utils_1.0.0        yaml_2.3.8              evaluate_0.24.0        
    ## INFO [2025-10-13 08:27:57]  [91] codetools_0.2-20        tibble_3.2.1            cli_3.6.3              
    ## INFO [2025-10-13 08:27:57]  [94] xtable_1.8-4            munsell_0.5.1           Rcpp_1.0.12            
    ## INFO [2025-10-13 08:27:57]  [97] gert_2.0.1              png_0.1-8               parallel_4.4.1         
    ## INFO [2025-10-13 08:27:57] [100] blob_1.2.4              doRNG_1.8.6             PERMANOVA_0.2.0        
    ## INFO [2025-10-13 08:27:57] [103] viridisLite_0.4.2       scales_1.3.0            openxlsx_4.2.8         
    ## INFO [2025-10-13 08:27:57] [106] purrr_1.0.2             crayon_1.5.3            clisymbols_1.2.0       
    ## INFO [2025-10-13 08:27:57] [109] GetoptLong_1.0.5        rlang_1.1.4             formatR_1.14
