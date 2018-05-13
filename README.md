# recountNNLSpaper

This repository contains the code used to derive the analyses presented in our paper.
links:
The scripts are presented in processing order.

01_annotation.R -> A R script that downloads the references used, and processes them for subsequent use.

02_indices.sh -> A shell script for downloading the genome as well as creating the alignment/quantification indices needed.

03_getPm.R -> A R script to estimate the probability matrices of each transcript to create observed feature counts.

04_get_bundles.R -> A R script to determine which genes need to be estimated together because they have overlapping mapping of reads.

05_simulate.R -> A R script that contains the code for the simulation studies of different methods using transcripts from protein-coding genes on chr1 and chr14.

06_CI.R -> A R script that contains the code for the simulation studies of the performance of our method's confidence intervals.

07_Geuvadis.sh -> A shell R script that details how we quantified the ERR188410 sample of the Geuvadis Consortium dataset using the methods covered in the paper.

08_simulate_rsem.R -> A R script that contains details how we evaluated the performance of the chosen models when ground truth is simulated using the RSEM estimates of sample ERR188410.

09_figures_tables.R -> A R script that creates the plots and tables presented (except the MA plots in the supplement).

10_graphics_MA.R -> A R script that produces the MA plots included in the supplement, as well as a short snippet to collect read lengths of samples in recount2.

11_run_all.R -> A R script that runs recountNNLS on all recount samples (minus TCGA and GTEX).

12_run_TCGA.R -> A R script that runs recountNNLS on the TCGA samples.

13_run_GTEX.R -> A R script that runs recountNNLS on the GTEX samples.

14_QCmetrics.R -> Short code snippet that we used checking how our fragment/read estimates correspond to overall mapped fragments/reads.

# Session info 
----------------------------------------------------------------------------------------------------------
 setting  value                                 
 version  R version 3.4.4 RC (2018-03-09 r74376)
 system   x86_64, linux-gnu                     
 ui       X11                                   
 language (EN)                                  
 collate  en_US.UTF-8                           
 tz       US/Eastern                            
 date     2018-05-10                            

Packages --------------------------------------------------------------------------------------------------------------
 package                     * version   date       source                            
 AnnotationDbi               * 1.40.0    2017-11-01 Bioconductor                      
 assertthat                    0.2.0     2017-04-11 CRAN (R 3.4.0)                    
 base                        * 3.4.4     2018-03-09 local                             
 Biobase                     * 2.38.0    2017-11-01 Bioconductor                      
 BiocGenerics                * 0.24.0    2017-11-01 Bioconductor                      
 BiocParallel                  1.12.0    2017-11-01 Bioconductor                      
 biomaRt                       2.34.2    2018-01-21 Bioconductor                      
 Biostrings                  * 2.46.0    2018-01-27 cran (@2.46.0)                    
 bit                           1.1-12    2014-04-09 CRAN (R 3.4.0)                    
 bit64                         0.9-7     2017-05-08 CRAN (R 3.4.1)                    
 bitops                        1.0-6     2013-08-17 CRAN (R 3.4.0)                    
 blob                          1.1.0     2017-06-17 CRAN (R 3.4.1)                    
 BSgenome                    * 1.46.0    2017-11-01 Bioconductor                      
 BSgenome.Hsapiens.UCSC.hg38 * 1.4.1     2017-05-26 Bioconductor                      
 compiler                      3.4.4     2018-03-09 local                             
 datasets                    * 3.4.4     2018-03-09 local                             
 DBI                           0.8       2018-03-02 CRAN (R 3.4.3)                    
 DelayedArray                * 0.4.1     2017-11-05 Bioconductor                      
 devtools                      1.13.5    2018-02-18 CRAN (R 3.4.3)                    
 digest                        0.6.15    2018-01-28 CRAN (R 3.4.3)                    
 GenomeInfoDb                * 1.14.0    2017-11-01 Bioconductor                      
 GenomeInfoDbData              1.0.0     2017-12-15 Bioconductor                      
 GenomicAlignments             1.14.1    2018-01-27 cran (@1.14.1)                    
 GenomicFeatures             * 1.30.3    2018-02-03 Bioconductor                      
 GenomicRanges               * 1.30.3    2018-02-27 Bioconductor                      
 graphics                    * 3.4.4     2018-03-09 local                             
 grDevices                   * 3.4.4     2018-03-09 local                             
 grid                          3.4.4     2018-03-09 local                             
 httr                          1.3.1     2017-08-20 CRAN (R 3.4.1)                    
 IRanges                     * 2.12.0    2017-11-01 Bioconductor                      
 lattice                       0.20-35   2017-03-25 CRAN (R 3.4.4)                    
 magrittr                      1.5       2014-11-22 CRAN (R 3.4.0)                    
 Matrix                        1.2-12    2017-11-30 CRAN (R 3.4.4)                    
 matrixStats                 * 0.53.1    2018-02-11 CRAN (R 3.4.3)                    
 memoise                       1.1.0     2017-04-21 CRAN (R 3.4.0)                    
 methods                     * 3.4.4     2018-03-09 local                             
 parallel                    * 3.4.4     2018-03-09 local                             
 pillar                        1.2.1     2018-02-27 CRAN (R 3.4.3)                    
 prettyunits                   1.0.2     2015-07-13 CRAN (R 3.4.0)                    
 progress                      1.1.2     2016-12-14 CRAN (R 3.4.0)                    
 R6                            2.2.2     2017-06-17 CRAN (R 3.4.1)                    
 Rcpp                          0.12.15   2018-01-20 cran (@0.12.15)                   
 RCurl                         1.95-4.10 2018-01-04 CRAN (R 3.4.3)                    
 recountNNLS                 * 0.99.7    2018-05-08 Github (jmf47/recountNNLS@9acad0b)
 recountNNLSdata             * 0.99.8    2018-02-21 Bioconductor                      
 rlang                         0.2.0     2018-02-20 CRAN (R 3.4.3)                    
 RMySQL                        0.10.14   2018-02-26 CRAN (R 3.4.3)                    
 Rsamtools                   * 1.30.0    2018-01-27 cran (@1.30.0)                    
 RSQLite                       2.0       2017-06-19 CRAN (R 3.4.1)                    
 rtracklayer                 * 1.38.3    2018-01-24 Bioconductor                      
 S4Vectors                   * 0.16.0    2017-11-01 Bioconductor                      
 stats                       * 3.4.4     2018-03-09 local                             
 stats4                      * 3.4.4     2018-03-09 local                             
 stringi                       1.1.6     2017-11-17 CRAN (R 3.4.2)                    
 stringr                     * 1.3.0     2018-02-19 cran (@1.3.0)                     
 SummarizedExperiment        * 1.8.1     2017-12-20 Bioconductor                      
 tibble                        1.4.2     2018-01-22 CRAN (R 3.4.3)                    
 tools                         3.4.4     2018-03-09 local                             
 utils                       * 3.4.4     2018-03-09 local                             
 withr                         2.1.1     2017-12-19 CRAN (R 3.4.3)                    
 XML                           3.98-1.10 2018-02-19 cran (@3.98-1.)                   
 XVector                     * 0.16.0    2018-05-10 Bioconductor                      
 zlibbioc                      1.24.0    2017-11-01 Bioconductor  