## Packages for  tami13

if(exists("course")){
    if(course == "EDO") level = 1
    if(course == "BE") level = 2
    }


pkg_bioc = NULL ## To accumulate Bioconductor packages
pkg_cran = NULL ## To accumulate CRAN packages


##  Basic packages installed using biocLite()
pkg_bioc_basico = c("Biobase","IRanges")
pkg_bioc = c(pkg_bioc,pkg_bioc_basico)

## Marginal differential expression
pkg_bioc_DE = NULL
if(level >= 1)
    pkg_bioc_DE = c("multtest", "genefilter","geneplotter","limma",
                    "Category","GOstats","qvalue","impute","oligo",
                    "golubEsets")
if(level >= 2)
    pkg_bioc_DE = c("siggenes")



pkg_bioc = c(pkg_bioc,pkg_bioc_DE)

## Data sets: microarrays
pkg_bioc_datasets_ma = c("ALLMLL","ALL")

pkg_bioc = c(pkg_bioc,pkg_bioc_datasets_ma)

## Affymetrix GeneChip
pkg_bioc_AffymetrixGeneChip = c("affy","simpleaffy","affycoretools","affyPLM",
                           "arrayQualityMetrics","gcrma")

pkg_bioc = c(pkg_bioc,pkg_bioc_AffymetrixGeneChip)

## Download data
pkg_bioc_DownloadData = c("GEOquery","ArrayExpress")

pkg_bioc = c(pkg_bioc,pkg_bioc_DownloadData)


## Gene Set Analysis  
pkg_bioc_GeneSetAnalysis = NULL
if(level >= 1)
    pkg_bioc_GeneSetAnalysis = c("GOstats","GSEABase","GSA",
                                 "EnrichmentBrowser")

if(level >= 3)
    pkg_bioc_GeneSetAnalysis = c("topGO","GSEAlm","PGSEA","GSVA",
                                 "goseq","globaltest")

pkg_bioc = c(pkg_bioc,pkg_bioc_GeneSetAnalysis)

## Network Analysis
pkg_bioc_NetworkAnalysis = NULL
if(level >= 3)
    pkg_bioc_NetworkAnalysis = c("DEGraph")

pkg_bioc = c(pkg_bioc,pkg_bioc_NetworkAnalysis)


## Annotation
pkg_bioc_annotation_generic = c("AnnotationDbi","biomaRt","rtracklayer",
                                "annotate","GO.db")

pkg_bioc_ChipDb = c("hgu133bcdf","hgu95av2.db","hgu133plus2.db","hgu133a.db",
                    "hgu133plus2probe",
                    "hgug4112a.db", ## Agilent 
                    "ath1121501.db", ## Affymetrix Arabidopsis ATH1 Genome Array
                    "pd.mogene.1.0.st.v1",
                    "ygs98.db")

pkg_bioc_OrgDb = NULL
if(level >= 1)
    pkg_bioc_OrgDb = c("org.Hs.eg.db")  ## Homo sapiens

if(level >= 2)
    pkg_bioc_OrgDb = c("org.Sc.sgd.db", ## yeast
                       "org.Mm.eg.db",  ## Mus Musculus
                       "org.Dm.eg.db",  ## Drosophila Melanogaster
                       "org.Rn.eg.db")


pkg_bioc_OrganismDb = NULL
if(level >= 1)
    pkg_bioc_OrganismDb = c("Homo.sapiens")
if(level >= 2)
    pkg_bioc_OrganismDb = c("Mus.musculus")

pkg_bioc_TxDb = NULL
if(level >= 2)
    pkg_bioc_TxDb = c("TxDb.Hsapiens.UCSC.hg19.knownGene",
                      "TxDb.Dmelanogaster.UCSC.dm3.ensGene",
                      "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")

pkg_bioc_BSgenome = NULL
if(level >= 2)
    pkg_bioc_BSgenome = c("BSgenome","BSgenome.Hsapiens.UCSC.hg19",
                          "BSgenome.Dmelanogaster.UCSC.dm3")

pkg_bioc_DataBase = NULL
if(level >= 2)
    pkg_bioc_DataBase = c("KEGGREST","geneLenDataBase")

pkg_bioc_annotation = c(pkg_bioc_annotation_generic,pkg_bioc_ChipDb,
                        pkg_bioc_OrgDb,pkg_bioc_OrganismDb,
                        pkg_bioc_TxDb,pkg_bioc_BSgenome,pkg_bioc_DataBase)


pkg_bioc = c(pkg_bioc,pkg_bioc_annotation)

## Reproducible research

pkg_bioc_reproducible_research = c("ReportingTools")
pkg_bioc = c(pkg_bioc,pkg_bioc_reproducible_research)

pkg_cran_reproducible_research = NULL
if(level >= 2)
    pkg_cran_reproducible_research =
        c("knitr","rmarkdown","reactable")

if(level  >= 3)
    pkg_cran_reproducible_research = c("HTMLUtils","R2HTML")

pkg_cran = c(pkg_cran,pkg_cran_reproducible_research)


pkg_cran = c(pkg_cran,pkg_cran_reproducible_research)

## Graphs
pkg_bioc_graphs = c("Rgraphviz")

pkg_bioc = c(pkg_bioc,pkg_bioc_graphs)

## microRNA
## pkg_bioc_microRNA = c("microRNA","mirna20.db","mirna20cdf","miRNApath",
##    "RmiR.Hs.miRNA","targetscan.Hs.eg.db","mirbase.db",
##                 "miRNAtap.db","RmiR.Hs.miRNA","RmiR.hsa")

## pkg_bioc = c(pkg_bioc,pkg_bioc_microRNA)

## Meta-analysis
##pkg_bioc_MetaAnalysis = c("MetaDE","GeneMeta","metaArray","RankProd","MAMA")


##pkg_bioc = c(pkg_bioc,pkg_bioc_MetaAnalysis) 


## RNASeq
pkg_bioc_RNASeq_basic = NULL


if(level >= 1)
pkg_bioc_RNASeq_basic = c("Biostrings","ShortRead","Rsamtools",
                          "GenomicRanges","GenomicAlignments",
                          "GenomicFeatures","SummarizedExperiment")

if(level >= 2)
pkg_bioc_RNASeq_basic = c("Gviz","EDASeq")

pkg_bioc_DatosRNASeq = NULL
if(level >= 2)
    pkg_bioc_DatosRNASeq = c("SRAdb","parathyroidSE","RTCGA")

pkg_bioc_RNASeqDE = NULL
if(level >= 1)
    pkg_bioc_RNASeqDE = c("edgeR")

if(level >= 2)
    pkg_bioc_RNASeqDE = c("DESeq2","parathyroidSE")

if(level >= 3)
    pkg_bioc_RNASeqDE = c("DESeq","NOISeq","goseq","cummeRbund")



pkg_bioc_scRNASeq = NULL
if(level >= 3)
    pkg_bioc_scRNASeq = c("Seurat","MiloR","cydar")

pkg_bioc_RNASeq = c(pkg_bioc_RNASeq_basic,pkg_bioc_DatosRNASeq,
                    pkg_bioc_RNASeqDE,pkg_bioc_scRNASeq)

pkg_bioc = c(pkg_bioc,pkg_bioc_RNASeq)

## Methylation 
pkg_bioc_methylation = NULL  
if(level >= 3)
pkg_bioc_methylation = c("minfi",
                         "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                         "IlluminaHumanMethylationEPICmanifest","missMethyl",
                         "minfiData","DMRcate")

pkg_bioc = c(pkg_bioc,pkg_bioc_methylation)

## Association studies
## To install SNPassoc previously
## apt install libfftw3-dev
## Bioconductor::BiocStyle

pkg_cran_association = pkg_bioc_association = NULL
if(level >= 3){
    pkg_cran_association = c("SNPassoc")
    pkg_bioc_association = c("BiocStyle")
}
pkg_cran = c(pkg_cran,pkg_cran_association)
pkg_bioc = c(pkg_bioc,pkg_bioc_association)


## Sequences
pkg_bioc_sequences = NULL
if(level >= 2)
    pkg_bioc_sequences = c("seqinr")

pkg_bioc = c(pkg_bioc,pkg_bioc_sequences)

## Proteomics
pkg_bioc_proteomics = NULL
if(level >= 3)
    pkg_bioc_proteomics = c("RforProteomics","mzR","mzID","MSnID","MSnbase",
                            "rpx","MLInterfaces","pRoloc","pRolocdata",
                            "MSGFplus","rols","hpar")

pkg_bioc = c(pkg_bioc,pkg_bioc_proteomics)

pkg_cran_proteomics = NULL
if(level >= 3)
    pkg_cran_proteomics = c("readMzXmlData")
pkg_cran = c(pkg_cran,pkg_cran_proteomics)


## Sample size and power
pkg_bioc_SampleSizePower = NULL
if(level >= 3)
    pkg_bioc_SampleSizePower = c("sizepower","OCplus","ssize",
                                 "RNASeqPower","RnaSeqSampleSize",
                                 "RnaSeqSampleSizeData","PROPER",
                                 "SSPA","GeneticsDesign","CSSP")

pkg_bioc = c(pkg_bioc,pkg_bioc_SampleSizePower)

## Tidy data tools
pkg_bioc_tidydata = NULL
if(level >= 2)
    pkg_bioc_tidydata = c("biobroom","ROC")

pkg_bioc = c(pkg_bioc,pkg_bioc_tidydata)


## R packages
pkg_cran_estadistica_basica = c("UsingR","nhstplot","DescTools")
pkg_cran = c(pkg_cran,pkg_cran_estadistica_basica)

pkg_cran_estadistica_multivariante = c("cluster")
pkg_cran = c(pkg_cran,pkg_cran_estadistica_multivariante)

pkg_cran_admin = NULL
if(level >= 2)
    pkg_cran_admin = c("devtools","roxygen2","pacman","pbapply")
pkg_cran = c(pkg_cran,pkg_cran_admin)

pkg_cran_utilidades = NULL
if(level >= 2)
    pkg_cran_utilidades = c("xtable","Rlab","R.matlab","Hmisc")
pkg_cran = c(pkg_cran,pkg_cran_utilidades)

pkg_cran_graphics = c("ggplot2","reshape","ggbio","ggfortify","ggdendro",
                      "factoextra","shinyFiles")
pkg_cran = c(pkg_cran,pkg_cran_graphics)


pkg_cran_varios = NULL
if(level >= 2)
    pkg_cran_varios = c("gplm","vcd","vcdExtra","spatstat","car","effects",
                        "gmodels", "fortunes")
pkg_cran = c(pkg_cran,pkg_cran_varios)


pkg_cran_bioinformatics = NULL
if(level >= 3)
    pkg_cran_bioinformatics = c("QuasiSeq","neat","kpmt")
pkg_cran = c(pkg_cran,pkg_cran_bioinformatics)

pkg_cran_tidydata = NULL
if(level >= 2) 
    pkg_cran_tidydata = c("broom")
pkg_cran = c(pkg_cran,pkg_cran_tidydata)


## Installation

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(pkg_bioc, version = "3.16")
install.packages(pkg_cran,repos = "https://cran.rediris.es/")


