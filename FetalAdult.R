#------phentypic information-----------------------
sample_run_phenotypic <- read.delim("C:/Users/Felix/Dropbox/DataScience/Projects/fetus_vs_adult_brains/gdc-sample-phenotypic-summalign_samples_analyzed.txt")

table2 <- sample_run_phenotypic[,c(1:11)]
library(pander)
pander(table2)

#-----Galaxy instance history-----------------------
library(GalaxyConnector) # load GalaxyConnector package 
# refer to https://github.com/scholtalbers/r-galaxy-connector for installation instructions

gx_init('b3eb4f05f812dc80752b4ea797753860', GALAXY_URL='https://usegalaxy.org/')

galaxy_hist <- gx_list_histories() # list all histories in the current Galaxy instance
galaxy_hist

#-----FastQ files------------------------------------
gx_init('b3eb4f05f812dc80752b4ea797753860', GALAXY_URL='https://usegalaxy.org/')
fastq_files <- gx_list_history_datasets() #assign all datasets to the fastq_files variable
fastq_files_new <- fastq_files[c(2:17), c(1:3)] # extracts datasets from 2 to 17, and print only 6 columns
pander(fastq_files_new)

#----Hisat BAM files----------------------------------
hisat2_job <- fastq_files[c(18:19, 26:31), 1:3] # extract metadata information of the hisat BAM files
pander(hisat2_job)

#----Alignment information and FastQC job----------------------------
library(plyr); library(dplyr) # load plyr and dplyr package to 
sample_run_phenotypic <- mutate (sample_run_phenotypic, pcntalign = ((align1 + alignp1)/input)*100)
pander(select(sample_run_phenotypic, run, ageg, input, align0, align1, alignp1, pcntalign))


fastqc_job <- fastq_files[33:48, 1:3]
pander(fastqc_job)

 <- gx_get(34, force = TRUE)
fastqc_541 <-gx_get(36, force = TRUE)
fastqc_566 <- gx_get(38, force = TRUE)
fastqc_567 <- gx_get(40, force = TRUE)

fastqc538 <- read.delim(fastqc_538, skip=12, nrows=55)
fastqc541 <- read.delim(fastqc_541, skip=12, nrows=55)
fastqc566 <- read.delim(fastqc_566, skip=12, nrows=55)
fastqc567 <- read.delim(fastqc_567, skip=12, nrows=55)

# download  raw fastqc test files of adult aligments
fastqc_534 <- gx_get(42, force = TRUE)
fastqc_535 <-gx_get(44, force = TRUE)
fastqc_536 <- gx_get(46, force = TRUE)
fastqc_561 <- gx_get(48, force = TRUE)

fastqc534 <- read.delim(fastqc_534, skip=12, nrows=55)
fastqc535 <- read.delim(fastqc_535, skip=12, nrows=55)
fastqc536 <- read.delim(fastqc_536, skip=12, nrows=55)
fastqc561 <- read.delim(fastqc_561, skip=12, nrows=55)

ageg <- factor(c(rep("fetus", 220), rep("adult", 220)))
ageg <- factor(ageg, order=TRUE, levels=c("fetus", "adult"))
fastqc_all <- data.frame(fastqc538$Mean, fastqc541$Mean, fastqc566$Mean, fastqc567$Mean, fastqc534$Mean, fastqc535$Mean,                     fastqc536$Mean, fastqc561$Mean)
names(fastqc_all) <-c("fqc538", "fqc541", "fqc566", "fqc567", "fqc534", "fqc535", "fqc536", "fqc561")

library(tidyr)
fastqc_all_tidy <- gather(fastqc_all, fastqc, pbsmeans, fqc538:fqc561)
ageg <- factor(c(rep("fetus", 220), rep("adult", 220)))
ageg <- factor(ageg, order=TRUE, levels=c("fetus", "adult"))
fastqc_all_final <- data.frame(fastqc_all_tidy, ageg)
fastqc_all_final$fastqc <- factor(fastqc_all_final$fastqc, order=TRUE, levels=c("fqc538", "fqc541", "fqc566", "fqc567", "fqc534", "fqc535", "fqc536", "fqc561"))

#-------Boxplot of FastQC----------------------------------
library(ggplot2)
fastqc_bxplot <- ggplot(fastqc_all_final, aes(x=fastqc, y=pbsmeans, fill=ageg)) + geom_boxplot()
fastqc_bxplot <- fastqc_bxplot + ggtitle("Per Base Sequence Quality by Hisat Runs and Age Groups") + 
            ylab("Average Quality Scores") + xlab("Hisat Alignment Runs")

fastqc_bxplot

fastqc_all_final %>% group_by (fastqc) %>% 
				summarise (mean_pbsq = mean(pbsmeans))

fastqc_age_bxplot <- ggplot(fastqc_all_final, aes(x=ageg, y=pbsmeans, fill=ageg)) + geom_boxplot()
fastqc_age_bxplot <- fastqc_age_bxplot + ggtitle("Per Base Sequence Quality by Age Group") + 
            ylab("Average Quality Scores") + xlab("Age Groups") + theme(legend.position="none")

fastqc_age_bxplot

fastqc_all_final %>% group_by (ageg) %>% 
				summarise (mean_pbsq = mean(pbsmeans))

#------feature counts--------------------------------
htseq_counts <- fastq_files[c(49, 51, 53, 55, 57, 59, 61, 63), c(1:3)]
pander(htseq_counts)

 <- read.delim(gx_get(49, force=TRUE))
htseq_541 <- read.delim(gx_get(51, force=TRUE))
htseq_566 <- read.delim(gx_get(53, force=TRUE))
htseq_567 <- read.delim(gx_get(55, force=TRUE))
htseq_534 <- read.delim(gx_get(57, force=TRUE))
htseq_535 <- read.delim(gx_get(59, force=TRUE))
htseq_536 <- read.delim(gx_get(61, force=TRUE))
htseq_561 <- read.delim(gx_get(63, force=TRUE))

gene_counts <- data.frame(htseq_538, htseq_541, htseq_566, htseq_567, htseq_534, htseq_535, htseq_536, htseq_561)
gene_counts_new <- gene_counts[,c(1:2,4,6,8,10,12, 14, 16)]
names(gene_counts_new) <- c("gene_id", "SRR1554538", "SRR1554541", "SRR1554566", "SRR1554567", "SRR1554534", "SRR1554535", "SRR1554536", "SRR1554561")

gene_counts_new1 <- select(gene_counts_new, gene_id, SRR1554538, SRR1554534, SRR1554541, SRR1554535, SRR1554566, SRR1554536, SRR1554567, SRR1554561)

pander(head(gene_counts_new1))

pander(tail(gene_counts_new1))

#------Exploratory analysis----------------------------
condition <- c("fetus", "adult", "fetus", "adult", "fetus", "adult", "fetus", "adult")
run_last3 <- c("538", "534", "541", "535", "566", "536", "567", "561")
col_data <- data.frame(condition, run_last3)
row.names(col_data) <- c("SRR1554538", "SRR1554534", "SRR1554541", "SRR1554535", "SRR1554566", "SRR1554536", "SRR1554567", "SRR1554561")
gene_counts_final <- gene_counts_new1[,2:9]
row.names(gene_counts_final) <- gene_counts_new1[,1]

library(DESeq2)
gene_dds <- DESeqDataSetFromMatrix(countData = gene_counts_final, 
            colData = col_data, design = ~ condition)

gene_dds_rltransf <- rlog(gene_dds, blind=FALSE)

dist_samples <- dist(t(assay(gene_dds_rltransf)))

gene_fit <- hclust(dist_samples, method="ward.D")

plot(gene_fit, hang=-1)

plotPCA(gene_dds_rltransf, intgroup = c("condition"))

#----Statistical Analysis------------------------------
gene_dd_dsq <- DESeq(gene_dds)

gene_exp_results <- results(gene_dd_dsq) #retrieve gene exp results
gene_exp_results_df <- as.data.frame(gene_exp_results) #results as data.frame

# sum of genes up- and down-regulated at > 1 or < -1 fold change, and padj < 0.5
sum(gene_exp_results_df$padj < 0.05 & gene_exp_results_df$log2FoldChange > 1, na.rm=TRUE)
sum(gene_exp_results_df$padj < 0.05 & gene_exp_results_df$log2FoldChange < -1, na.rm=TRUE)

## generate volcano plot
with(gene_exp_results_df, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

## add color to points 
with(subset(gene_exp_results_df, padj<.05 ), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(gene_exp_results_df, abs(log2FoldChange)>1), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(gene_exp_results_df, padj<.05 & abs(log2FoldChange)>1), 
     points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

library (dplyr) # version 0.5.0
gene_exp_results_df <- data.frame(gene_counts_new[,1], gene_exp_results_df)
gene_exp_results_df %>% filter (log2FoldChange > 1) %>% arrange(padj) %>% head # top up-reg genes
gene_exp_results_df %>% filter (log2FoldChange < -1) %>% arrange(padj) %>% head # top down-reg genes

#-----Gene set analysis----------------------------------
(AnnotationHub) # version 2.4.2
# assign annotations records to the ah
ah <- AnnotationHub()

### perform query to ID files corresonding to fetal brain samples and H3K4me3 
fetal_brain <- ah[["AH44720"]]; # fetus dataset
adult_brain <- ah[["AH43565"]]; # adult dataset
liver <- ah[["AH44167"]] ## liver dataset

# evaluate the widths between peaks for each sample
summary(width(fetal_brain)); summary(width(adult_brain)); summary(width(liver))

# generate a VennDiagrama to visualize overlapping peaks between 2 of the samples
## or all three samples at once
library(ChIPpeakAnno) # version 3.6.5
ol_f_a_l <- findOverlapsOfPeaks(fetal_brain, adult_brain, liver)

makeVennDiagram(ol_f_a_l)

#---Gene set analysis part II----------------------------
ref_seq <- query(ah, "RefSeq")
ref_seq_hg19 <- ref_seq[ref_seq$genome == "hg19" & ref_seq$title== "RefSeq Genes"]
ref_seq_hg19 <- ref_seq_hg19[[1]] ## download the information

## calculate percent of promoters with H3K4me3 peaks
promoters_ref_seq_hg19 <- promoters(ref_seq_hg19) 

fetal_brain_ovlps <- findOverlaps(promoters_ref_seq_hg19, 
                                  fetal_brain) # overlaps of fetal brain promoters
adult_brain_ovlps <- findOverlaps(promoters_ref_seq_hg19, 
                                  adult_brain) ## find overlaps of adult brain promoters
liver_ovlps <- findOverlaps(promoters_ref_seq_hg19, liver) # find overlaps of liver promoters

### calculate the percentages of promoters with H3K4me3 peaks
fetal_prcnt_pm_pk <- length(unique(subjectHits(fetal_brain_ovlps))) / length(promoters_ref_seq_hg19)
adult_prcnt_pm_pk <- length(unique(subjectHits(adult_brain_ovlps))) / length(promoters_ref_seq_hg19)
liver_prcnt_pm_pk <- length(unique(subjectHits(liver_ovlps))) / length(promoters_ref_seq_hg19)

prcnt_pm_pk_df <- data.frame(fetal_prcnt_pm_pk, adult_prcnt_pm_pk, liver_prcnt_pm_pk)
names(prcnt_pm_pk_df) <-c("fetal", "adult", "liver")

library(tidyr) # version 0.6.0
prcnt_pm_pk_df_tidy <- gather(prcnt_pm_pk_df, sample, H3K4me3_pcnt_pks, fetal:liver)
pander(prcnt_pm_pk_df_tidy)