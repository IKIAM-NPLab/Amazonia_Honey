Spectral deconvolution of stingless bee honeys HS-GC-MS metabolomics
================
Angiely Camacho, José Abata, Jefferson Pastuna
2024-04-12

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#before-to-start" id="toc-before-to-start">Before to start</a>
- <a href="#erah-package-workflow" id="toc-erah-package-workflow">eRah
  package workflow</a>
  - <a href="#compound-deconvolution"
    id="toc-compound-deconvolution">Compound deconvolution</a>
  - <a href="#alignment" id="toc-alignment">Alignment</a>
  - <a href="#missing-compound-recovery"
    id="toc-missing-compound-recovery">Missing compound recovery</a>
- <a href="#identification" id="toc-identification">Identification</a>

# Introduction

The present document aims to record the HS-GC-EI (Q)MS spectral
deconvolution procedure of the volatile profile of honeys of 3 native
bee species (*Melipona fasciculata*, *Melipona fuscopilosa*, and
*Tetragonisca angustula*) from 3 different localities in the Chontapunta
parish, Tena, Ecuadorian Amazonia region.

The workflow used is taken from the paper [eRah: A Computational Tool
Integrating Spectral Deconvolution and Alignment with Quantification and
Identification of Metabolites in GC/MS-Based
Metabolomics](https://doi.org/10.1021/acs.analchem.6b02927). It offers a
wide variety of functions to automatically detect and deconvolve the
spectra of the compounds appearing in GC–MS chromatograms. The “*eRah*”
R package includes automatic peak detection, peak alignment, feature
extraction, and compound identification.

# Before to start

The “*eRah*” R package accepts raw data files (netCDF or mzXML) obtained
in GC–q/MS, GC-TOF/MS, and GC-qTOF/MS (using nominal mass) equipment.

In this case, the data was obtained with a GC-2030 gas chromatograph
coupled to a GCMS-QP2020 NX quadrupole mass spectrometer operating in
electron ionization (GC-EI-q/MS). The raw MS data (.qgd) were converted
to (netCDF) format using the proprietary software of the instrument GCMS
Postrun Analysis 4.53SP1. The data was organized in an experiment folder
named Data_to_eRah. This contains three class-folders called ’Blank’,
’Quality_Control’, and ’Samples’, each containing the sample files for
that class.

# eRah package workflow

The “*eRah*” R package workflow is based on five steps: (i) data
pre-processing, (ii) spectral deconvolution, (iii) spectral alignment,
(iv) retrieval of missing compounds, and (v) compound identification.

To use the “*eRah*” R package, you need to install the “*eRah*” package
and call the library() function.

``` r
# eRah package installation
#install.packages('erah')
# eRah library call
library(erah)
```

To proceed, you must delete the unwanted files in a specific path and
create a directory with the wanted files.

``` r
# Delete all file that are not in folders
unlink('Data/Data_to_eRah/*')
# Data folder path
createdt('Data/Data_to_eRah/')
```

Load and process the necessary data from the CSV files, which contain
the relevant information to proceed with the study about the volatilome
of honey from different Amazonian stingless bees. The MetaboSet object
was created to store the sample metadata and process.

**NOTE.** The instrumental and phenotype files created by eRah were
relocated to a new directory. In both files, the semicolon delimiter has
been changed to a comma delimiter. The complete directory path of the
chromatograms has been added to the instrumental file.

``` r
# Loading (*.mzXML) chromatograms name (instrumental file)
#instrumental <- read.csv('Data/Metadata_to_eRah/Metadata_inst_mzXML.csv')
# If (*.mzXML) did not work
# Loading (*.CDF) chromatograms
instrumental <- read.csv('Data/Metadata_to_eRah/HS_GCMS_Data_inst_CDF.csv')
# Loading metadata of the chromatograms (phenotype file)
phenotype <- read.csv('Data/Metadata_to_eRah/HS_GCMS_Data_pheno.csv')
# Merge of metadata information with chromatograms
ex <- newExp(instrumental = instrumental,
             phenotype = phenotype,
             info = 'Amazonia honey')
```

## Compound deconvolution

It is a mathematical and computational technique used to recover an
original signal from a complex pool of information.

Set up the parameters with the minimum height and width of the peak, the
noise threshold, and the mz signals compounds wanted to exclude.

``` r
ex.dec.par <- setDecPar(min.peak.width = 2,
                        min.peak.height = 450,
                        noise.threshold = 45,
                        avoid.processing.mz = c(30:69,73:75,147:149))
```

To run parallel processing, the package “future” must be used.

**NOTE.** Adjusts the number of “workers” with the equipped PC cores.

``` r
plan(future::multisession,
     workers = 16)
```

We proceed to the deconvolution of compounds using the parameters
specified in “dec_par”. The results obtained will be saved in
(dec_data).

``` r
ex <- deconvolveComp(ex,
                     ex.dec.par)
```

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Blank/4_Blanco.CDF ... Processing 1 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/15_QC1_4.CDF ... Processing 2 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/21_Qc2_2.CDF ... Processing 3 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/28_Qc3_3.CDF ... Processing 4 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/34_Qc1_2.CDF ... Processing 5 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/41_Qc2_4.CDF ... Processing 6 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/44_Qc2_5.CDF ... Processing 7 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/45_Qc3_2.CDF ... Processing 8 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/7_Qc1_5.CDF ... Processing 9 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Quality_Control/8_Qc2_3.CDF ... Processing 10 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/10_RIG003_1.CDF ... Processing 11 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/11_Sample_FIED001_2.CDF ... Processing 12 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/12_Sample_RGY001_3.CDF ... Processing 13 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/13_Sample_RIG003_3.CDF ... Processing 14 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/14_Sample_RYTA006_2.CDF ... Processing 15 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/16_Sample_RYMG001_3.CDF ... Processing 16 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/17_Sample_RGY001_1.CDF ... Processing 17 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/18_Sample_FIED001_1.CDF ... Processing 18 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/19_Sample_RYMG001_2.CDF ... Processing 19 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/20_Sample_ER001_1.CDF ... Processing 20 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/22_Sample_ATA001_2.CDF ... Processing 21 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/23_RIG003_2.CDF ... Processing 22 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/25_Sample_RYTA006_3.CDF ... Processing 23 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/26_Sample_RYTA006_1.CDF ... Processing 24 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/27_Sample_RIG005_3.CDF ... Processing 25 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/29_Sample_RGB004_1.CDF ... Processing 26 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/30_Sample_RGB004_3.CDF ... Processing 27 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/31_Sample_FIED001_3.CDF ... Processing 28 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/32_Sample_ATA001_3.CDF ... Processing 29 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/33_Sample_ATA001_1.CDF ... Processing 30 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/35_Sample_RGY001_2.CDF ... Processing 31 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/36_Sample_RGB004_2.CDF ... Processing 32 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/37_Sample_RIG005_1.CDF ... Processing 33 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/38_Sample_ER001_2.CDF ... Processing 34 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/40_Sample_RYMG001_1.CDF ... Processing 35 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/42_Sample_RIG005_2.CDF ... Processing 36 / 37

    ## 
    ##  Deconvolving compounds from Data/Data_to_eRah/Samples/43_Sample_ER001_3.CDF ... Processing 37 / 37

    ## Compounds deconvolved

## Alignment

Correcting for retention time variation of eluting compounds,
facilitating relative quantification and comparison of compounds between
samples, is done by grouping compounds within these limits of retention
time distance and spectral similarity.

``` r
# Alignment parameters
ex.al.par <- setAlPar(min.spectra.cor = 0.80,
                      max.time.dist = 3,
                      mz.range = 70:550)
# Alignment
ex <- alignComp(ex,
                alParameters = ex.al.par)
```

## Missing compound recovery

Missing Compound Recovery is a function that aims to identify and
recover information about compounds that may have been missed or
undetected in the initial analysis due to various reasons, such as
background noise, low signal intensity, or overlapping peaks.

To recover missing compounds, use the function “recMissComp” and save
the result in “ex.”

``` r
ex <- recMissComp(ex,
                  min.samples = 3)
```

    ## 
    ##  Updating alignment table... 
    ## Model fitted!

Exporting alignment feature list.

``` r
# Extracting alignment feature list
feat_list <- alignList(ex,
                       by.area = FALSE)
# Exporting alignment feature list
#write.csv(feat_list,
#          file = "Result/eRah_Result/erah_Export_2Notame.csv")
```

# Identification

Metabolite identification was by comparing all the spectra found against
a reference database. The “*eRah*” R package default database was used
for metabolite identification. Metabolite identification was improved by
exporting all the spectra found (.msp) in the “*eRah*” R package to the
NIST MS Search 2.4 software.

``` r
# Identification
peak_iden <- identifyComp(ex,
                          id.database = mslib,
                          mz.range = NULL,
                          n.putative = 1)
```

    ## Constructing matrix database... 
    ## Comparing spectra... 
    ## Done!

Exporting spectra to NIST MS Search software for identification with the
NIST-20 library.

``` r
#export2MSP(peak_iden,
#           store.path = "Result/eRah_Result",
#           alg.version = 2)
```
