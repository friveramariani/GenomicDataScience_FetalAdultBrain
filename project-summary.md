#Title: RNA-seq re-analysis to demonstrate students the RNA-seq data analysis workflow

## Summary
In this genomic data science project, which served as a demonstration for students in an upper division microbiology course (MCB3023 at Miami Dade College), the differential gene expression between fetal and adult brains was evaluated. Transcriptome sequencing data (known as RNA-seq) for this genomic data science project, which has been previously analyzed and published [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4281298/), was sequenced on an Illumina platform from human post-mortem brains. The raw sequencing data related to each brain were uploaded from a **[public database](http://www.ebi.ac.uk/ena/data/view/PRJNA245228)** into the **Galaxy Project server** (https://usegalaxy.org), while the meta-data from each sample was obtained manually from the 3rd link listed below. The data was aligned, quality controls performed for each alignment, and gene count levels of expression determined. Furthermore, exploratory analysis and fitting statistical models were elaborated to identify patterns of gene expression based on phenotypic variables of the RNA-seq samples. From the statistical model, 18 up-regulated and 18 down-regulated genes were subsetted and assigned to students for them to search the function of these genes. In-class discussions were carried out to elaborate the significance of the genes functions and the importance of up or down-regulating the respective genes. 

### Sites to access the 1) original article, 2)  RNA-seq data, 3) phenotype meta-data, and 4) how to get the data 4

1. Find in this link the related publication: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4281298/

2. Find in this link the RNA-seq data: http://www.ebi.ac.uk/ena/data/view/PRJNA245228

3. Find in this link the phenotype meta-data for the samples: http://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA245228

### Steps performed
1. Download the RNA-seq data and phenotype meta-data
2. RNA-seq sequences algined through **The Galaxy Project, https://usegalaxy.org**)
3. QC of alignments
4. Calculate gene expressions based on gene count levels
5. Perform exploratory analysis to refine model to be built
6. Built model to ID genes differentially expressed between fetal and adult brains
7. Identify possible biological patterns
8. Assign genes to students
9. Document the whole process through a brief report




