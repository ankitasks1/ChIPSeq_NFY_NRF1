print("Identification of Allele Specific Peaks")
print("Factor in Analysis :  NFY")
print("Mouse Genome : version mm10")
print("Mouse cell line : BJ1")
print("Raw data processing pipeline : /media/ankitv/Archivio1/2021/chip-seq/BJ1_NFY")

library(DESeq2)
library(EDASeq)
library(metaseqR)
library(splitstackshape)
library(ggplot2)
library(plyr)
library(BaalChIP)
setwd("/media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/enrichment_R1")
#But I used macs1 output becuase it has slightly higher number of peaks
cp /media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/bulk/NFY/rep1/PeakCall_R1_MACS1/NFY_rep1_peaks.bed ./
cp /media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/bulk/NFY/rep2/PeakCall_R1_MACS1/NFY_rep2_peaks.bed ./

#Count Peaks
wc -l NFY_rep1_peaks.bed 
wc -l NFY_rep2_peaks.bed

#Index the bams which can be used in visualization
samtools index ./../bulk/Input/Input_R1_uniq.rmdup.bam 
samtools index ./../bulk/NFY/rep1/NFY_rep1_R1_uniq.rmdup.bam 
samtools index ./../bulk/NFY/rep2/NFY_rep2_R1_uniq.rmdup.bam 
samtools index ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam
samtools index ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam
samtools index ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam




#Bulk peak intersection (using only -a and -b option will aloow to intersect peaks in such a way that only shared region between two peaks will be printed)
bedtools intersect -a NFY_rep1_peaks.bed -b NFY_rep2_peaks.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"NR"\t"$5}' > NFY_intersection_rep1_2_peaks.bed
#I selected No -u and let shared  intervals present and given NUMBER to each peak ID
#Coverage calculation, 3 possible ways
bedtools multicov -bams ./../bulk/Input/Input_R1_uniq.rmdup.bam -bed NFY_intersection_rep1_2_peaks.bed > Input_R1_NFY_int_peaks.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../bulk/Input/Input_R1_uniq.rmdup.bam  > Input_NFY_intersection_rep1_2_peaks.bed
#and featurecounts based
#For Bulk already contain chr in chromosme files

#But allel specific data does not contain chr so chr need to be added before peak calling and featurecounts
#First convert all picard dedup bam to bed (NOTE: When I convert picard output bed '.rmdup.' is used while manualy dedup file has '_rmdup.'. So this is the difference and can be seen in mnauscript)
#eg. Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed (picard)
#Input_rep1_R1_uniq.sortedByReadname.genome.sort_JF1_rmdup.bed

#Run BamtoBed in respective folders
bedtools bamtobed -i Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
bedtools bamtobed -i Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed


bedtools bamtobed -i NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
bedtools bamtobed -i NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed

bedtools bamtobed -i NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
bedtools bamtobed -i NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed


awk '{print $4"\t"$1"\t"$2"\t"$3"\t""."}' NFY_intersection_rep1_2_peaks.bed > NFY_intersection_rep1_2_peaks.saf
#/home/ankitv/tools_av/subread-2.0.0-source/bin/featureCounts -a NFY_intersection_rep1_2_peaks.saf -F SAF -o BJ1_NFY_featurecounts_R1.txt ./../bulk/Input/Input_R1_uniq.rmdup.bam  ./../bulk/NFY/rep1/NFY_rep1_R1_uniq.rmdup.bam  ./../bulk/NFY/rep2/NFY_rep2_R1_uniq.rmdup.bam  ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam  ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam

NFY_intersection_rep1_2_peaks <- read.table("NFY_intersection_rep1_2_peaks.bed", header = F)
colnames(NFY_intersection_rep1_2_peaks) <- c("chr","start","end","peakID","score")


#So I used bedtools coverage based read assignment
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../bulk/Input/Input_R1_uniq.rmdup.bam > Input_R1_uniq.rmdup.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../bulk/NFY/rep1/NFY_rep1_R1_uniq.rmdup.bam > NFY_rep1_R1_uniq.rmdup.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../bulk/NFY/rep2/NFY_rep2_R1_uniq.rmdup.bam > NFY_rep2_R1_uniq.rmdup.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/NFY/rep1/NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b ./../allele_sp/NFY/rep2/NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed

paste Input_R1_uniq.rmdup.NFY_cov.bed NFY_rep1_R1_uniq.rmdup.NFY_cov.bed NFY_rep2_R1_uniq.rmdup.NFY_cov.bed Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed NFY_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NFY_cov.bed NFY_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NFY_cov.bed > BJ1_NFY_covdata_R1.txt

#testfilefeat <- read.table("BJ1_NFY_featurecounts_R1_re.txt")
#colnames(testfilefeat) <- c("peakID","chr","start","end","strand","length","Input","NFY_1","NFY_2","Input_B6","Input_JF1","NFY_B6_1","NFY_JF1_1","NFY_B6_2","NFY_JF1_2")


#So I used bedtools coverage based read assignment
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1 <- read.table("BJ1_NFY_covdata_R1.txt", header = F)
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1)
dim(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1)
#Test if paste parallely
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(seq(4,length(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1),9))])
tail(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(seq(4,length(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1),9))])
#So the data is pasted parallely
#Now prepare coverage file
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1 <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(1:5,seq(6,length(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1),9))]
colnames(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1) <- c("chr","start","end","peakID","score","Input","NFY_1","NFY_2","Input_B6","Input_JF1","NFY_B6_1","NFY_JF1_1","NFY_B6_2","NFY_JF1_2")
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
tail(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
#match with head(testfilefeat) tail((testfilefeat))

head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
tail(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
dim(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
rownames(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1) <- paste0(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1$chr,"%",
                                                                   NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1$start,"%",
                                                                   NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1$end,"%",
                                                                   NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1$peakID,"%",
                                                                   NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1$score)
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1 <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,c(6:14)]

head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
colSums(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr["peaks"] <- rownames(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- cSplit(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr, "peaks", "%") #library(splitstackshape)
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
dim(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr[,c(10:14, 1:9)]
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
write.table(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr,"NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Use Deseq2 for calculating scaling factor from bulk
#Input normalization not required as the peaks are called on background of input
head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
Bulk_NFY_countdata <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,2:3]
head(Bulk_NFY_countdata)


bulkNFYcoldata <- read.table("bulkNFYcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bulkNFYcoldata) <- colnames(Bulk_NFY_countdata)
head(bulkNFYcoldata)
bulkNFYcoldata <- bulkNFYcoldata[,c("condition","replicate")]
bulkNFYcoldata$condition <- factor(bulkNFYcoldata$condition)
bulkNFYcoldata$replicate <- factor(bulkNFYcoldata$replicate)
all(rownames(bulkNFYcoldata) == colnames(Bulk_NFY_countdata)) #should print TRUE
#Low counts filter (filter only rowsums for NFY allele specific)
#NFYkeep <- rowSums(Bulk_NFY_countdata) >= 10
#Bulk_NFY_countdatafilt <- Bulk_NFY_countdata[NFYkeep,]
#dim(Bulk_NFY_countdatafilt)

#Filter using proportion function of NOISeq:

Replicates <- factor(bulkNFYcoldata$replicate)
dim(Bulk_NFY_countdata)
Bulk_NFY_countdatafilt <- NOISeq::filtered.data(dataset=Bulk_NFY_countdata, factor=Replicates, method=3, norm=FALSE)
head(Bulk_NFY_countdatafilt)
plotRLE(as.matrix(Bulk_NFY_countdatafilt), outline=F)
write.table(Bulk_NFY_countdatafilt, "Bulk_NFY_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)


bulkNFYdds <- DESeqDataSetFromMatrix(countData =Bulk_NFY_countdatafilt, colData = bulkNFYcoldata, design = ~ replicate)
bulkNFYdds
#No filtering
bulkNFYddscounts <- counts(bulkNFYdds, normalized=FALSE)
head(bulkNFYddscounts)
dim(bulkNFYddscounts)
colSums(bulkNFYddscounts)
#View filtered count matrix: View(counts(bulkNFYdds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
bulkNFYddsNorm <- estimateSizeFactors(bulkNFYdds)
sizeFactors(bulkNFYddsNorm)

head(NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
Alsp_NFY_countdata <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,6:9]
head(Alsp_NFY_countdata)
gcoldata <- read.table("NFYgcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(gcoldata)
gcoldata[,1]
rownames(gcoldata)=gcoldata[,1]
rownames(gcoldata)
colnames(gcoldata)
gcoldata = gcoldata[,-1]
gcoldata <- data.frame(gcoldata)
head(gcoldata)
gcoldata <- gcoldata[,c("gcondition","greplicate","gallele")]
rownames(gcoldata) <- colnames(Alsp_NFY_countdata)
gcoldata$gcondition <- factor(gcoldata$gcondition)
gcoldata$greplicate <- factor(gcoldata$greplicate)
gcoldata$gallele <- factor(gcoldata$gallele)

all(rownames(gcoldata) == colnames(Alsp_NFY_countdata)) #should print TRUE

#Low counts filter (filter only rowsums for NFY allele specific)
#alNFYkeep <- rowSums(Alsp_NFY_countdata[,6:9]) >= 10
#Alsp_NFY_countdatafilt <- Alsp_NFY_countdata[alNFYkeep,]
ASReplicates <- factor(gcoldata$gallele)
dim(Alsp_NFY_countdata)
Alsp_NFY_countdatafilt <- NOISeq::filtered.data(dataset=Alsp_NFY_countdata, factor=ASReplicates, method=3, norm=FALSE)
head(Alsp_NFY_countdatafilt)
plotRLE(as.matrix(Alsp_NFY_countdatafilt), outline=F)
write.table(Alsp_NFY_countdatafilt, "Alsp_NFY_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

dim(Alsp_NFY_countdatafilt)
AlspNFYdds <- DESeqDataSetFromMatrix(countData =Alsp_NFY_countdatafilt, colData = gcoldata, design = ~ gallele)
AlspNFYdds
#No filtering
AlspNFYddscounts <- counts(AlspNFYdds, normalized=FALSE)
head(AlspNFYddscounts)
dim(AlspNFYddscounts)
colSums(AlspNFYddscounts)

#Normalize by DESeq2 size factor, set size factoe using sizeFactors(bulkNFYddsNorm)
AlspNFYddsNorm <- AlspNFYdds
sizeFactors(AlspNFYddsNorm) <- c(0.942809, 0.942809, 1.060660, 1.060660)
sizeFactors(AlspNFYddsNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
AlspNFYddsNormcounts <- counts(AlspNFYddsNorm, normalized=TRUE)
head(AlspNFYddsNormcounts)
write.table(AlspNFYddsNormcounts, "AlspNFYddsNormcounts.txt", sep="\t", quote=F, append = F)
boxplot(AlspNFYddsNormcounts, ylim=c(0,100))
AlspNFYddsNormcounts <- data.frame(AlspNFYddsNormcounts)
AlspNFYddsNormcounts["matAllelicRatio1"] <- (AlspNFYddsNormcounts$NFY_JF1_1 / (AlspNFYddsNormcounts$NFY_B6_1 + AlspNFYddsNormcounts$NFY_JF1_1))
AlspNFYddsNormcounts["matAllelicRatio2"] <- (AlspNFYddsNormcounts$NFY_JF1_2 / (AlspNFYddsNormcounts$NFY_B6_2 + AlspNFYddsNormcounts$NFY_JF1_2))
head(AlspNFYddsNormcounts)
dim(AlspNFYddsNormcounts)

AlspNFYddsNormcounts_chr <- AlspNFYddsNormcounts
AlspNFYddsNormcounts_chr["peaks"] <- rownames(AlspNFYddsNormcounts_chr)
head(AlspNFYddsNormcounts_chr)
AlspNFYddsNormcounts_chr <- cSplit(AlspNFYddsNormcounts_chr, "peaks", "%") #library(splitstackshape)
head(AlspNFYddsNormcounts_chr)
AlspNFYddsNormcounts_chr <- AlspNFYddsNormcounts_chr[,c(7:11, 1:6)]
head(AlspNFYddsNormcounts_chr)
colnames(AlspNFYddsNormcounts_chr) <- c("chr", "start", "end", "peakID", "score","NFY_B6_1","NFY_JF1_1","NFY_B6_2","NFY_JF1_2","matAllelicRatio1","matAllelicRatio2")
write.table(AlspNFYddsNormcounts_chr,"AlspNFYddsNormcounts_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


##Predict Monoallelic and Biallelic Peaks  --------------###
#------------ I used Proportion test -------------#
#Combined p.vales -2Sumation(log(Pi))
#ChiSquare df = 2*length(p-values), so if 2 p.values 2 *2
#lower.tail =FALSE and lower.tail =TRUE, 1-p
#Prop test and tag chromsome positions
#AlspNFYddsNormcounts_chr <- read.table("AlspNFYddsNormcounts_chr.txt", header = F)

head(AlspNFYddsNormcounts_chr)
dim(AlspNFYddsNormcounts_chr)
AlspNFYddsNormcounts_chr <- data.frame(AlspNFYddsNormcounts_chr)
rownames(AlspNFYddsNormcounts_chr) <- paste(AlspNFYddsNormcounts_chr$chr,AlspNFYddsNormcounts_chr$start,AlspNFYddsNormcounts_chr$end,AlspNFYddsNormcounts_chr$peakID,AlspNFYddsNormcounts_chr$score, sep = "%")
AlspNFYddsNormcounts_chr_prop <- AlspNFYddsNormcounts_chr[,6:11]
head(AlspNFYddsNormcounts_chr_prop)
#Take sum of alleles
AlspNFYddsNormcounts_chr_prop["NFY_alavg_1"] <- AlspNFYddsNormcounts_chr_prop$NFY_B6_1 + AlspNFYddsNormcounts_chr_prop$NFY_JF1_1
AlspNFYddsNormcounts_chr_prop["NFY_alavg_2"] <- AlspNFYddsNormcounts_chr_prop$NFY_B6_2 + AlspNFYddsNormcounts_chr_prop$NFY_JF1_2

#Take allelic ratio
AlspNFYddsNormcounts_chr_prop["mat_allele_ratio"] <- (AlspNFYddsNormcounts_chr_prop$matAllelicRatio1 + AlspNFYddsNormcounts_chr_prop$matAllelicRatio2)/2


dim(AlspNFYddsNormcounts_chr_prop) #4770    9
AlspNFYddsNormcounts_chr_propid <- AlspNFYddsNormcounts_chr_prop
AlspNFYddsNormcounts_chr_propid["id"] <- rownames(AlspNFYddsNormcounts_chr_prop)
head(AlspNFYddsNormcounts_chr_propid)
writexl::write_xlsx(AlspNFYddsNormcounts_chr_propid, "AlspNFYddsNormcounts_chr_prop.xlsx")

dim(AlspNFYddsNormcounts_chr_prop)#4770    9
head(AlspNFYddsNormcounts_chr_prop,1)
summary(AlspNFYddsNormcounts_chr_prop)
AlspNFYddsNormcounts_chr_prop <- data.frame(AlspNFYddsNormcounts_chr_prop)
#NA detected as some in some cases both allele were 0, So That need to be removed
AlspNFYddsNormcounts_chr_prop.filt <- na.omit(AlspNFYddsNormcounts_chr_prop) 
head(AlspNFYddsNormcounts_chr_prop.filt)
dim(AlspNFYddsNormcounts_chr_prop.filt)
summary(AlspNFYddsNormcounts_chr_prop.filt)
NFY_alavg_1.prop1 <- Map(prop.test,x =AlspNFYddsNormcounts_chr_prop.filt$NFY_B6_1, n= AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_1, p=0.5)

AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg_1.pvalue"] <- data.frame(capture.output(for (i in seq_along(NFY_alavg_1.prop1)){
  cat(NFY_alavg_1.prop1[[i]]$p.value, "\n")
}))

NFY_alavg_2.prop1 <- Map(prop.test,x =AlspNFYddsNormcounts_chr_prop.filt$NFY_B6_2, n= AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_2, p=0.5)

AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg_2.pvalue"] <- data.frame(capture.output(for (i in seq_along(NFY_alavg_2.prop1)){
  cat(NFY_alavg_2.prop1[[i]]$p.value, "\n")
}))


#Combine p-values calculated from prop.test for rep1 and rep2
AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg.comb.testpvalue"] <- -2 *(log(as.numeric(as.character(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_1.pvalue))) + log(as.numeric(as.character(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_2.pvalue))))
#To deal with -2 * (log(0)+log(0)), I used fisher.method function which has option to zero.sub to value, I set it 2.2e-16, this will replace Inf with 2.2e-16, 
#eg. -2*(log(2.2e-16)+log(2.2e-16)) = 144.2116
class(AlspNFYddsNormcounts_chr_prop.filt)
AlspNFYddsNormcounts_chr_prop.filt <- data.frame(AlspNFYddsNormcounts_chr_prop.filt)
AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg.comb.fisherpvalue"] <- fisher.method(data.frame(as.numeric(as.character(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_1.pvalue)), as.numeric(as.character(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg_2.pvalue))), method = c("fisher"), p.corr ="none", zero.sub = 2.2e-16, na.rm = FALSE, mc.cores=NULL)
#Warning will pop up Warning message:
#In `[<-.data.frame`(`*tmp*`, "NFY_alavg.comb.fisherpvalue", value = list( :provided 4 variables to replace 1 variables
#Export the sheet and manually checked the column NFY_alavg.comb.testpvalue and NFY_alavg.comb.fisherpvalue should match except Inf, so manual and fisher.method is ok
head(AlspNFYddsNormcounts_chr_prop.filt)
AlspNFYddsNormcounts_chr_prop.filt_check <- AlspNFYddsNormcounts_chr_prop.filt
AlspNFYddsNormcounts_chr_prop.filt_check["ID"] <- rownames(AlspNFYddsNormcounts_chr_prop.filt_check)
writexl::write_xlsx(AlspNFYddsNormcounts_chr_prop.filt_check, "AlspNFYddsNormcounts_chr_prop.filt_check.xlsx")

head(AlspNFYddsNormcounts_chr_prop.filt)
tail(AlspNFYddsNormcounts_chr_prop.filt)
#Chisquare test for p-values#pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)#p vector of probabilities, n	 number of observations, degrees of freedom (non-negative, but can be non-integer)
AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg.comb.pvalue.pchisq"] <- pchisq(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg.comb.fisherpvalue,4, lower.tail=FALSE)
AlspNFYddsNormcounts_chr_prop.filt["NFY_alavg.comb.fdr"] <- p.adjust(AlspNFYddsNormcounts_chr_prop.filt$NFY_alavg.comb.pvalue.pchisq,method="BH")
head(AlspNFYddsNormcounts_chr_prop.filt)

#Introduce peak info
AlspNFYddsNormcounts_chr_prop.filt["id"] <- data.frame(rownames(AlspNFYddsNormcounts_chr_prop.filt))
head(AlspNFYddsNormcounts_chr_prop.filt)
AlspNFYddsNormcounts_chr_prop.filt.sep <- cSplit(AlspNFYddsNormcounts_chr_prop.filt, "id", "%")
dim(AlspNFYddsNormcounts_chr_prop.filt.sep)
head(AlspNFYddsNormcounts_chr_prop.filt.sep)
AlspNFYddsNormcounts_chr_prop.filt.sep <- AlspNFYddsNormcounts_chr_prop.filt.sep[,c(16:20,1:15)]
colnames(AlspNFYddsNormcounts_chr_prop.filt.sep) <- c("chr","start","end","peakID","score",colnames(AlspNFYddsNormcounts_chr_prop.filt[,1:15]))
#JF1_monoallelic
AlspNFYddsNormcounts_JF1_mono <- AlspNFYddsNormcounts_chr_prop.filt.sep[which(AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.80 &
                                                                                AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.80),]
head(AlspNFYddsNormcounts_JF1_mono)
dim(AlspNFYddsNormcounts_JF1_mono)


#Monoalleically expressed to JF1 with FDR < 0.05
AlspNFYddsNormcounts_JF1_monofilt <- AlspNFYddsNormcounts_JF1_mono[which(AlspNFYddsNormcounts_JF1_mono$NFY_alavg.comb.fdr<0.05),]
head(AlspNFYddsNormcounts_JF1_monofilt)
dim(AlspNFYddsNormcounts_JF1_monofilt)
write.table(AlspNFYddsNormcounts_JF1_monofilt,"AlspNFYddsNormcounts_JF1_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNFYddsNormcounts_JF1_monofilt_qual <- AlspNFYddsNormcounts_JF1_monofilt
head(AlspNFYddsNormcounts_JF1_monofilt_qual)
AlspNFYddsNormcounts_JF1_monofilt_qual199 <- AlspNFYddsNormcounts_JF1_monofilt_qual[which(AlspNFYddsNormcounts_JF1_monofilt_qual$score >= 199),]
dim(AlspNFYddsNormcounts_JF1_monofilt_qual199)
head(AlspNFYddsNormcounts_JF1_monofilt_qual199)
write.table(AlspNFYddsNormcounts_JF1_monofilt_qual199,"AlspNFYddsNormcounts_JF1_monofilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#B6_monoallelic
AlspNFYddsNormcounts_B6_mono <- AlspNFYddsNormcounts_chr_prop.filt.sep[which(AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.20 &
                                                                               AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.20),]
head(AlspNFYddsNormcounts_B6_mono)
dim(AlspNFYddsNormcounts_B6_mono)

#Monoalleically expressed to B6 with FDR < 0.05
AlspNFYddsNormcounts_B6_monofilt <- AlspNFYddsNormcounts_B6_mono[which(AlspNFYddsNormcounts_B6_mono$NFY_alavg.comb.fdr<0.05),]
head(AlspNFYddsNormcounts_B6_monofilt)
dim(AlspNFYddsNormcounts_B6_monofilt)
write.table(AlspNFYddsNormcounts_B6_monofilt,"AlspNFYddsNormcounts_B6_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNFYddsNormcounts_B6_monofilt_qual <- AlspNFYddsNormcounts_B6_monofilt
head(AlspNFYddsNormcounts_B6_monofilt_qual)
AlspNFYddsNormcounts_B6_monofilt_qual199 <- AlspNFYddsNormcounts_B6_monofilt_qual[which(AlspNFYddsNormcounts_B6_monofilt_qual$score >= 199),]
dim(AlspNFYddsNormcounts_B6_monofilt_qual199)
head(AlspNFYddsNormcounts_B6_monofilt_qual199)
write.table(AlspNFYddsNormcounts_B6_monofilt_qual199,"AlspNFYddsNormcounts_B6_monofilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#Biallelic
#Biallelicalleleic
AlspNFYddsNormcounts_Biallelic <- AlspNFYddsNormcounts_chr_prop.filt.sep[which(AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.60 &
                                                                                 AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.60 &
                                                                                 AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.40 &
                                                                                 AlspNFYddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.40),]
head(AlspNFYddsNormcounts_Biallelic)
dim(AlspNFYddsNormcounts_Biallelic)

#Bialleically expressed  with FDR > 0.05
AlspNFYddsNormcounts_Biallelicfilt <- AlspNFYddsNormcounts_Biallelic[which(AlspNFYddsNormcounts_Biallelic$NFY_alavg.comb.fdr>0.05),]
head(AlspNFYddsNormcounts_Biallelicfilt)
dim(AlspNFYddsNormcounts_Biallelicfilt)
write.table(AlspNFYddsNormcounts_Biallelicfilt,"AlspNFYddsNormcounts_Biallelicfilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNFYddsNormcounts_Biallelicfilt_qual <- AlspNFYddsNormcounts_Biallelicfilt
head(AlspNFYddsNormcounts_Biallelicfilt_qual)
AlspNFYddsNormcounts_Biallelicfilt_qual199 <- AlspNFYddsNormcounts_Biallelicfilt_qual[which(AlspNFYddsNormcounts_Biallelicfilt_qual$score >= 199),]
dim(AlspNFYddsNormcounts_Biallelicfilt_qual199)
head(AlspNFYddsNormcounts_Biallelicfilt_qual199)
write.table(AlspNFYddsNormcounts_Biallelicfilt_qual199,"AlspNFYddsNormcounts_Biallelicfilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

##No need to overlap on allele specific called peaks
##bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt.txt -b ./../merge_peakanalysis/NFY_alsp_uniqB6_bulk.bed | sort -k4,4 -u > AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt
##bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt.txt -b ./../merge_peakanalysis/NFY_alsp_uniqJF1_bulk.bed | sort -k4,4 -u > AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt
##bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt.txt -b ./../merge_peakanalysis/NFY_alsp_sharedB6JF1_bulk.bed | sort -k4,4 -u > AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt

cp AlspNFYddsNormcounts_B6_monofilt.txt AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt
cp AlspNFYddsNormcounts_JF1_monofilt.txt AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt
cp AlspNFYddsNormcounts_Biallelicfilt.txt AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt

##Overlap with ICR mm10 (for checking monoallelecity) 
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt


#Perform similar analysis on Input
#Purpose: Heterozygosity map using Input of Bulk and Allele specific
### ----------------Heterozygous Peaks in Input---------------------- ###

Bulk_Input_countdata <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,c(1,1)]
head(Bulk_Input_countdata)
Alsp_Input_countdata <- NFY_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,4:5]
head(Alsp_Input_countdata)

bulkInputcoldata <- read.table("bulkInputcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bulkInputcoldata) <- colnames(Bulk_Input_countdata)
head(bulkInputcoldata)
bulkInputcoldata <- bulkInputcoldata[,c("condition","replicate")]
bulkInputcoldata$condition <- factor(bulkInputcoldata$condition)
bulkInputcoldata$replicate <- factor(bulkInputcoldata$replicate)
all(rownames(bulkInputcoldata) == colnames(Bulk_Input_countdata)) #should print TRUE
#Low counts filter (filter only rowsums for Input allele specific)
#Inputkeep <- rowSums(Bulk_Input_countdata) >= 10
#Bulk_Input_countdatafilt <- Bulk_Input_countdata[Inputkeep,]
#dim(Bulk_Input_countdatafilt)

#Filter using proportion function of NOISeq:

InReplicates <- factor(bulkInputcoldata$replicate)
dim(Bulk_Input_countdata)
Bulk_Input_countdatafilt <- NOISeq::filtered.data(dataset=Bulk_Input_countdata, factor=InReplicates, method=3, norm=FALSE)
head(Bulk_Input_countdatafilt)
plotRLE(as.matrix(Bulk_Input_countdatafilt), outline=F)
#write.table(Bulk_Input_countdatafilt, "Bulk_Input_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)


bulkInputdds <- DESeqDataSetFromMatrix(countData =Bulk_Input_countdatafilt, colData = bulkInputcoldata, design = ~ replicate)
bulkInputdds
#No filtering
bulkInputddscounts <- counts(bulkInputdds, normalized=FALSE)
head(bulkInputddscounts)
dim(bulkInputddscounts)
colSums(bulkInputddscounts)
#View filtered count matrix: View(counts(bulkInputdds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
bulkInputddsNorm <- estimateSizeFactors(bulkInputdds)
sizeFactors(bulkInputddsNorm)

head(Alsp_Input_countdata)
Inputgcoldata <- read.table("Inputgcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(Inputgcoldata)
Inputgcoldata[,1]
rownames(Inputgcoldata)=Inputgcoldata[,1]
rownames(Inputgcoldata)
colnames(Inputgcoldata)
Inputgcoldata = Inputgcoldata[,-1]
Inputgcoldata <- data.frame(Inputgcoldata)
head(Inputgcoldata)
Inputgcoldata <- Inputgcoldata[,c("gcondition","greplicate","gallele")]
rownames(Inputgcoldata) <- colnames(Alsp_Input_countdata)
Inputgcoldata$gcondition <- factor(Inputgcoldata$gcondition)
Inputgcoldata$greplicate <- factor(Inputgcoldata$greplicate)
Inputgcoldata$gallele <- factor(Inputgcoldata$gallele)

all(rownames(Inputgcoldata) == colnames(Alsp_Input_countdata)) #should print TRUE

#Low counts filter (filter only rowsums for Input allele specific)
#alInputkeep <- rowSums(Alsp_Input_countdata[,6:9]) >= 10
#Alsp_Input_countdatafilt <- Alsp_Input_countdata[alInputkeep,]
ASInReplicates <- factor(Inputgcoldata$gallele)
dim(Alsp_Input_countdata)
Alsp_Input_countdatafilt <- NOISeq::filtered.data(dataset=Alsp_Input_countdata, factor=ASInReplicates, method=3, norm=FALSE)
head(Alsp_Input_countdatafilt)
plotRLE(as.matrix(Alsp_Input_countdatafilt), outline=F)
#write.table(Alsp_Input_countdatafilt, "Alsp_Input_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

dim(Alsp_Input_countdatafilt)
AlspInputdds <- DESeqDataSetFromMatrix(countData =Alsp_Input_countdatafilt, colData = Inputgcoldata, design = ~ gallele)
AlspInputdds
#No filtering
AlspInputddscounts <- counts(AlspInputdds, normalized=FALSE)
head(AlspInputddscounts)
dim(AlspInputddscounts)
colSums(AlspInputddscounts)

#Normalize by DESeq2 size factor, set size factoe using sizeFactors(bulkInputddsNorm)
AlspInputddsNorm <- AlspInputdds
sizeFactors(AlspInputddsNorm) <- c(1, 1)
sizeFactors(AlspInputddsNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
AlspInputddsNormcounts <- counts(AlspInputddsNorm, normalized=TRUE)
head(AlspInputddsNormcounts)
write.table(AlspInputddsNormcounts, "AlspInputddsNormcounts.txt", sep="\t", quote=F, append = F)
boxplot(AlspInputddsNormcounts, ylim=c(0,100))
AlspInputddsNormcounts <- data.frame(AlspInputddsNormcounts)
AlspInputddsNormcounts["matAllelicRatio"] <- (AlspInputddsNormcounts$Input_JF1 / (AlspInputddsNormcounts$Input_B6 + AlspInputddsNormcounts$Input_JF1))
head(AlspInputddsNormcounts)
dim(AlspInputddsNormcounts)

AlspInputddsNormcounts_chr <- AlspInputddsNormcounts
AlspInputddsNormcounts_chr["peaks"] <- rownames(AlspInputddsNormcounts_chr)
head(AlspInputddsNormcounts_chr)
AlspInputddsNormcounts_chr <- cSplit(AlspInputddsNormcounts_chr, "peaks", "%") #library(splitstackshape)
head(AlspInputddsNormcounts_chr)
AlspInputddsNormcounts_chr <- AlspInputddsNormcounts_chr[,c(4:8, 1:3)]
head(AlspInputddsNormcounts_chr)
colnames(AlspInputddsNormcounts_chr) <- c("chr", "start", "end", "peakID", "score","Input_B6","Input_JF1","matAllelicRatio")
write.table(AlspInputddsNormcounts_chr,"AlspInputddsNormcounts_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#--------------------   Identify Monoallelic and Biallelic peaks in Input ----------------------#

#------------  Proportion test -------------#
#Combined p.vales -2Sumation(log(Pi))
#ChiSquare df = 2*length(p-values), so if 2 p.values 2 *2
#lower.tail =FALSE and lower.tail =TRUE, 1-p
#Prop test and tag chromsome positions
#AlspInputddsNormcounts_chr <- read.table("AlspInputddsNormcounts_chr.txt", header = F)

head(AlspInputddsNormcounts_chr)
dim(AlspInputddsNormcounts_chr)
AlspInputddsNormcounts_chr <- data.frame(AlspInputddsNormcounts_chr)
rownames(AlspInputddsNormcounts_chr) <- paste(AlspInputddsNormcounts_chr$chr,AlspInputddsNormcounts_chr$start,AlspInputddsNormcounts_chr$end,AlspInputddsNormcounts_chr$peakID,AlspInputddsNormcounts_chr$score, sep = "%")
AlspInputddsNormcounts_chr_prop <- AlspInputddsNormcounts_chr[,6:8]
head(AlspInputddsNormcounts_chr_prop)
#Take sum of alleles
AlspInputddsNormcounts_chr_prop["Input_alavg"] <- AlspInputddsNormcounts_chr_prop$Input_B6 + AlspInputddsNormcounts_chr_prop$Input_JF1

#Take allelic ratio

dim(AlspInputddsNormcounts_chr_prop) #4827    4
AlspInputddsNormcounts_chr_propid <- AlspInputddsNormcounts_chr_prop
AlspInputddsNormcounts_chr_propid["id"] <- rownames(AlspInputddsNormcounts_chr_prop)
head(AlspInputddsNormcounts_chr_propid)
#writexl::write_xlsx(AlspInputddsNormcounts_chr_propid, "AlspInputddsNormcounts_chr_prop.xlsx")

dim(AlspInputddsNormcounts_chr_prop)#4827    4
head(AlspInputddsNormcounts_chr_prop,1)
summary(AlspInputddsNormcounts_chr_prop)
AlspInputddsNormcounts_chr_prop <- data.frame(AlspInputddsNormcounts_chr_prop)
#NA detected as some in some cases both allele were 0, So That need to be removed
AlspInputddsNormcounts_chr_prop.filt <- na.omit(AlspInputddsNormcounts_chr_prop) 
head(AlspInputddsNormcounts_chr_prop.filt)
dim(AlspInputddsNormcounts_chr_prop.filt)
summary(AlspInputddsNormcounts_chr_prop.filt)
Input_alavg_1.prop1 <- Map(prop.test,x =AlspInputddsNormcounts_chr_prop.filt$Input_B6, n= AlspInputddsNormcounts_chr_prop.filt$Input_alavg, p=0.5)

AlspInputddsNormcounts_chr_prop.filt["Input_alavg_1.pvalue"] <- data.frame(capture.output(for (i in seq_along(Input_alavg_1.prop1)){
  cat(Input_alavg_1.prop1[[i]]$p.value, "\n")
}))

AlspInputddsNormcounts_chr_prop.filt["Input_alavg.comb.fdr"] <- p.adjust(AlspInputddsNormcounts_chr_prop.filt$Input_alavg_1.pvalue,method="BH")
head(AlspInputddsNormcounts_chr_prop.filt)

#Introduce gene name
AlspInputddsNormcounts_chr_prop.filt["id"] <- data.frame(rownames(AlspInputddsNormcounts_chr_prop.filt))
head(AlspInputddsNormcounts_chr_prop.filt)
AlspInputddsNormcounts_chr_prop.filt.sep <- cSplit(AlspInputddsNormcounts_chr_prop.filt, "id", "%")
dim(AlspInputddsNormcounts_chr_prop.filt.sep)
head(AlspInputddsNormcounts_chr_prop.filt.sep)
AlspInputddsNormcounts_chr_prop.filt.sep <- AlspInputddsNormcounts_chr_prop.filt.sep[,c(7:11,1:6)]
colnames(AlspInputddsNormcounts_chr_prop.filt.sep) <- c("chr","start","end","peakID","score",colnames(AlspInputddsNormcounts_chr_prop.filt[,1:6]))
#JF1_monoallelic
AlspInputddsNormcounts_JF1_mono <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio >= 0.80),]
head(AlspInputddsNormcounts_JF1_mono)
dim(AlspInputddsNormcounts_JF1_mono)


#Monoalleically expressed to JF1 with FDR < 0.05
AlspInputddsNormcounts_JF1_monofilt <- AlspInputddsNormcounts_JF1_mono[which(AlspInputddsNormcounts_JF1_mono$Input_alavg.comb.fdr<0.05),]
head(AlspInputddsNormcounts_JF1_monofilt)
dim(AlspInputddsNormcounts_JF1_monofilt)
write.table(AlspInputddsNormcounts_JF1_monofilt,"AlspInputddsNormcounts_JF1_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#B6_monoallelic
AlspInputddsNormcounts_B6_mono <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio <= 0.20),]
head(AlspInputddsNormcounts_B6_mono)
dim(AlspInputddsNormcounts_B6_mono)

#Monoalleically expressed to B6 with FDR < 0.05
AlspInputddsNormcounts_B6_monofilt <- AlspInputddsNormcounts_B6_mono[which(AlspInputddsNormcounts_B6_mono$Input_alavg.comb.fdr<0.05),]
head(AlspInputddsNormcounts_B6_monofilt)
dim(AlspInputddsNormcounts_B6_monofilt)
write.table(AlspInputddsNormcounts_B6_monofilt,"AlspInputddsNormcounts_B6_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)
#Biallelic
#Biallelicalleleic
AlspInputddsNormcounts_Biallelic <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio <= 0.60 &
                                                                                     AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio >= 0.40),]
head(AlspInputddsNormcounts_Biallelic)
dim(AlspInputddsNormcounts_Biallelic)

#Bialleically expressed  with FDR > 0.05
AlspInputddsNormcounts_Biallelicfilt <- AlspInputddsNormcounts_Biallelic[which(AlspInputddsNormcounts_Biallelic$Input_alavg.comb.fdr>0.05),]
head(AlspInputddsNormcounts_Biallelicfilt)
dim(AlspInputddsNormcounts_Biallelicfilt)
write.table(AlspInputddsNormcounts_Biallelicfilt,"AlspInputddsNormcounts_Biallelicfilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)



##Identify Motif covered by SNP 
#Monoallelic B6
#Create a motif overlapped file Black 6
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/mouse/motif[CCAAT].mm10.bed CCAAT_motif.mm10.sort.bed
#CCAAT and JF1-B6 SNP
bedtools intersect -wa -wb -a CCAAT_motif.mm10.sort.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > CCAAT_jf1v2SNPmm10.bed

#Motif and SNP overlap peaks
#CCAAT
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b CCAAT_motif.mm10.sort.bed > NFY_B6_specific_peaks_CCAAT.bed

#SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_B6_specific_peaks_jf1v2.bed

#Note CCAAT_jf1v2SNPmm10.bed is file intersected with SNP (jf1v2mm10)
#CCAAT+SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b CCAAT_jf1v2SNPmm10.bed > NFY_B6_specific_peaks_CCAAT_jf1v2.bed

#CCAAT, SNP
bedtools intersect -wa -wb -a NFY_B6_specific_peaks_CCAAT.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_B6_specific_peaks_CCAAT_then_jf1v2.bed

#NO Nfy
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b CCAAT_motif.mm10.sort.bed -v > Nomotif_NFY_B6_specific_peaks.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b CCAAT_jf1v2SNPmm10.bed -v > NomotifSNP_NFY_B6_specific_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NFY_B6_specific_peaks.bed

#Monoallelic JF1
#Create a motif overlapped file JF1
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/mouse/motif[CCAAT].mm10.jf1.bed CCAAT_motif.mm10.jf1.sort.bed
#CCAAT and JF1-JF1 SNP
bedtools intersect -wa -wb -a CCAAT_motif.mm10.jf1.sort.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > CCAAT_jf1v2SNPmm10.jf1.bed

#Motif and SNP overlap peaks
#CCAAT
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b CCAAT_motif.mm10.jf1.sort.bed > NFY_JF1_specific_peaks_CCAAT.bed

#SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_JF1_specific_peaks_jf1v2.bed

#Note CCAAT_jf1v2SNPmm10.jf1.bed is file intersected with SNP (jf1v2mm10.jf1)
#CCAAT+SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b CCAAT_jf1v2SNPmm10.jf1.bed > NFY_JF1_specific_peaks_CCAAT_jf1v2.bed

#CCAAT, SNP
bedtools intersect -wa -wb -a NFY_JF1_specific_peaks_CCAAT.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_JF1_specific_peaks_CCAAT_then_jf1v2.bed

#NO Nfy
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b CCAAT_motif.mm10.jf1.sort.bed -v > Nomotif_NFY_JF1_specific_peaks.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b CCAAT_jf1v2SNPmm10.jf1.bed -v > NomotifSNP_NFY_JF1_specific_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NFY_JF1_specific_peaks.bed

#Indel
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndel.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel.txt

#CCAAT and JF1-B6 Indels
bedtools intersect -wa -wb -a CCAAT_motif.mm10.sort.bed -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed > CCAAT_jf1vindelSNPmm10.bed
bedtools intersect -wa -wb -a CCAAT_motif.mm10.jf1.sort.bed -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed > CCAAT_jf1vindelSNPmm10.jf1.bed

#CCAAT+Indels
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b CCAAT_jf1vindelSNPmm10.bed > NFY_B6_specific_peaks_CCAAT_jf1vindel.bed

#CCAAT+Indels
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b CCAAT_jf1vindelSNPmm10.jf1.bed > NFY_JF1_specific_peaks_CCAAT_jf1vindel.bed

#4321634 /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed
#32233 CCAAT_jf1vindelSNPmm10.bed
#32045 CCAAT_jf1vindelSNPmm10.jf1.bed


#------------------------- Identify Heterozygous SNPs  ------------------------#

filepath <- "/media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/enrichment_R1"
awk '{print $1,$2-100,$3+100,$4,$5}' OFS='\t' NFY_rep1_summits.bed > NFY_rep1_summits_span100.bed
awk '{print $1,$2-100,$3+100,$4,$5}' OFS='\t' NFY_rep2_summits.bed > NFY_rep2_summits_span100.bed

#Share span extended summits in rep1 and rep2
bedtools intersect -a NFY_rep1_summits_span100.bed -b NFY_rep2_summits_span100.bed  > NFY_intersection_rep1_2_summits_span100.bed 

#Get SNPs overlapped summit +/- 100bp
bedtools intersect -wa -wb -a NFY_intersection_rep1_2_summits_span100.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed | awk '{print $6,$7,$8,$9}' OFS='_' | sort -k1,1 -u | awk -F'_' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'  > jf1v2_Snp_GRm10_NFYsummitsp100.bed
awk '{print NR"\t"$1"\t"$2"\t"$4"\t"$5}' jf1v2_Snp_GRm10_NFYsummitsp100.bed > jf1v2_hetSNP.txt
vim jf1v2_hetSNP.txt
ID	CHROM	POS	REF	ALT
jf1v2hets <- c("BJ1"="jf1v2_hetSNP.txt")
head(read.delim(file.path(filepath,"jf1v2_hetSNP.txt")))
#Copy summit span to other name
cp NFY_intersection_rep1_2_summits_span100.bed BJ1_NFY.bed
mysamplesheet <- file.path(filepath,"baalchipmysamplesheet.tsv")

#mysamplesheet <- read.delim(file.path(filepath,"mysamplesheet.tsv"))
mysamplesheet

BJ1res <- BaalChIP(samplesheet=mysamplesheet, hets=jf1v2hets)
BJ1res
BaalChIP.get(BJ1res, what="samples")
#run alleleCounts
BJ1res <- alleleCounts(BJ1res, min_base_quality=10, min_mapq=15, verbose=FALSE)

#run QC filter
#res <- QCfilter(res, 
#                RegionsToFilter=list("blacklist"=blacklist_mm10), 
#                verbose=FALSE)
BJ1res <- mergePerGroup(BJ1res)
BJ1res <- filter1allele(BJ1res)
BJ1res <- getASB(BJ1res, Iter=5000, conf_level=0.95, cores = 8, 
                 RMcorrection = FALSE, 
                 RAFcorrection=FALSE)
Sys.time()
BJ1result <- BaalChIP.report(BJ1res)
head(BJ1result[["BJ1"]])
#show ASB SNPs
BJ1resultInput <- data.frame(BJ1result[["BJ1"]])
head(BJ1resultInput)
dim(BJ1resultInput)
BJ1resultInput$REF.counts <- BJ1resultInput$REF.counts/2
BJ1resultInput$ALT.counts <- BJ1resultInput$ALT.counts/2
BJ1resultInput$Total.counts <- BJ1resultInput$Total.counts/2
BJ1resultInput <- BJ1resultInput[which(BJ1resultInput$REF.counts >= 2),]
BJ1resultInput <- BJ1resultInput[which(BJ1resultInput$ALT.counts >= 2),]

write.table(BJ1resultInput, "BJ1resultInputBi.txt", sep = "\t", append = F, quote = F, row.names = F)

BJ1resultInputASBSNP <- BJ1resultInput[BJ1resultInput$isASB==TRUE,]
head(BJ1resultInputASBSNP)
dim(BJ1resultInputASBSNP)
write.table(BJ1resultInputASBSNP, "BJ1resultInputASBSNP.txt", sep = "\t", append = F, quote = F, row.names = F)


BJ1resultInputBi <- read.table("BJ1resultInputBi.txt", header = T, stringsAsFactors = F)
head(BJ1resultInputBi)
dim(BJ1resultInputBi)
BJ1resultInputBi_chr <- BJ1resultInputBi[,c(2,3,3:15)]
head(BJ1resultInputBi_chr)
write.table(BJ1resultInputBi_chr,  "BJ1resultInputBi_chr.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


#Heterozygous SNPs file: BJ1resultInputBi_chr.txt
#Heterozygouse/Biallleic region in Input  file: AlspInputddsNormcounts_Biallelicfilt.txt
#Add biases peak in input information to matrix 

bedtools intersect -wa -wb -a AlspInputddsNormcounts_Biallelicfilt.txt -b BJ1resultInputBi_chr.txt > AlspInputddsNormcounts_Biallelicfilt_BJ1resultInputBi.txt

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b BJ1resultInputBi_chr.txt > AlspNFYddsNormcounts_B6_BJ1InputBi.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b BJ1resultInputBi_chr.txt > AlspNFYddsNormcounts_JF1_BJ1InputBi.txt

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNFYddsNormcounts_B6_InputBiallelicfilt.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNFYddsNormcounts_JF1_InputBiallelicfilt.txt


bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b AlspInputddsNormcounts_Biallelicfilt_BJ1resultInputBi.txt > AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b AlspInputddsNormcounts_Biallelicfilt_BJ1resultInputBi.txt > AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi.txt


#------------ The further characteristics of peaks are related to regions +/- 100 from the middle of peaks
#Peaks at +/- 100bp from center of peaks
AlspNFYddsNormcounts_B6_monofilt_uniqB6 <- read.table("AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100 <- AlspNFYddsNormcounts_B6_monofilt_uniqB6
head(AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100)
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100["diff"] <- round((AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$V3-AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$V2)/2)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100["chr"] <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$V1
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100["start"] <- (AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$V2 + AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$diff)-100
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100["end"] <- (AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$V3 - AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100$diff)+100
head(AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100 <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100[,c(22:24,4,5)]
write.table(AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100, "AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)

AlspNFYddsNormcounts_JF1_monofilt_uniqJF1 <- read.table("AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100 <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1
head(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100["diff"] <- round((AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V3-AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V2)/2)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100["chr"] <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V1
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100["start"] <- (AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V2 + AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$diff)-100
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100["end"] <- (AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V3 - AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100$diff)+100
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100 <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100[,c(22:24,4,5)]
write.table(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100, "AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)

AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1 <- read.table("AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100 <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1
head(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["diff"] <- round((AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V3-AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V2)/2)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["chr"] <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V1
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["start"] <- (AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V2 + AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$diff)-100
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["end"] <- (AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V3 - AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$diff)+100
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100 <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100[,c(22:24,4,5)]
write.table(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100, "AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)




#----------- Get peaks status on CpG Islands  ----------------#
CGIs_mm10 <- read.table("CGIs_mm10.txt", header = F, stringsAsFactors = F)
head(CGIs_mm10)
CGIs_mm10_chr <- CGIs_mm10[,-c(1,7:12)]
head(CGIs_mm10_chr)
write.table(CGIs_mm10_chr, "CGIs_mm10_chr.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

#Overlap with NFY extended peak center
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NFY_B6_specific_peaks_CGIs_mm10_chr.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chr.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr.txt


#---------------------------  Correlation with Cofactors  -------------------------------#
#Get Cofactors data
cat jf1v2_Snp+chr.GRCm38.mm10.bed  /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38_ucsc_mm10.bed | sort -k1,1 -k2,2n > jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed
#Only files are stored
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Sp1/Sp1_motif.mm10.bed  -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > SP1_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Sp1/Sp1_peaks.bed  -b SP1_motif_SNP_Indels_mm10.bed > SP1_peaks_motif_SNP_Indels_mm10.bed


bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Yy1/YY1_motif.mm10.bed  -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > YY1_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Yy1/DMSO_YY1_rep1_peaks.bed  -b YY1_motif_SNP_Indels_mm10.bed > YY1_peaks_motif_SNP_Indels_mm10.bed


#NRF1 the files are taken from newway analysis
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NRF1/rep1/PeakCall_R1_MACS1/NRF1_rep1_peaks.bed ./
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NRF1/rep2/PeakCall_R1_MACS1/NRF1_rep2_peaks.bed ./
  
  
bedtools intersect -a NRF1_rep1_peaks.bed -b NRF1_rep2_peaks.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"NR"\t"$5}' > NRF1_intersection_rep1_2_peaks.bed
bedtools intersect -wa -wb -a NRF1_motif.mm10.sort.bed -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > NRF1_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a NRF1_intersection_rep1_2_peaks.bed -b NRF1_motif_SNP_Indels_mm10.bed > NRF1_peaks_motif_SNP_Indels_mm10.bed



#Overlap NFY peaks with feature

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NFY_B6_specific_peaks_SP1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NFY_JF1_specific_peaks_SP1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels.bed

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NFY_B6_specific_peaks_YY1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NFY_JF1_specific_peaks_YY1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels.bed

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b NRF1_peaks_motif_SNP_Indels_mm10.bed > NFY_B6_specific_peaks_NRF1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b NRF1_peaks_motif_SNP_Indels_mm10.bed > NFY_JF1_specific_peaks_NRF1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b NRF1_peaks_motif_SNP_Indels_mm10.bed > NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indels.bed

#------------------------ Public datasets --------------------------------#
#All files from NFYA, NFYB, NFYC were combined to get merged files and lifovered to mm10
#/media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NFYs/lift_mm10/NFY_mm10.bed

#Overlap NFY 
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt -b /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NFYs/lift_mm10/NFY_mm10.bed > NFY_B6_specific_peaks_pub_NFY.bed

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt -b /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NFYs/lift_mm10/NFY_mm10.bed > NFY_JF1_specific_peaks_pub_NFY.bed

bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b  /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NFYs/lift_mm10/NFY_mm10.bed  > NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY.bed


#---------------------  Countings --------------------#
#Intersections count
awk '{print $1"_"$2"_"$3"_"$4}' CCAAT_jf1v2SNPmm10.bed | sort -k1,1 -u | wc -l
awk '{print $1"_"$2"_"$3"_"$4}' CCAAT_jf1v2SNPmm10.jf1.bed | sort -k1,1 -u | wc -l

#NFY_B6
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt_qual199.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt | wc -l
sort -k4,4 -u NFY_B6_specific_peaks_CCAAT.bed | wc -l
sort -k4,4 -u NFY_B6_specific_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NFY_B6_specific_peaks_CCAAT_jf1v2.bed | wc -l
sort -k4,4 -u NFY_B6_specific_peaks_CCAAT_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NFY_B6_specific_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NFY_B6_specific_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NFY_B6_specific_peaks.bed | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_BJ1InputBi.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_InputBiallelicfilt.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_B6_monofilt_uniqB6_midspan100.txt | wc -l

#NFY_JF1
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt_qual199.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt | wc -l
sort -k4,4 -u NFY_JF1_specific_peaks_CCAAT.bed | wc -l
sort -k4,4 -u NFY_JF1_specific_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NFY_JF1_specific_peaks_CCAAT_jf1v2.bed | wc -l
sort -k4,4 -u NFY_JF1_specific_peaks_CCAAT_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NFY_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NFY_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NFY_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_BJ1InputBi.txt  | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_InputBiallelicfilt.txt  | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi.txt  | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt  | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt | wc -l

#--------------------------------    Quality filter ------------------------------#
awk '{if($5>=199) print $0}' AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt > AlspNFYddsNormcounts_B6_monofilt_uniqB6.qual199.txt
awk '{if($5>=199) print $0}' AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt > AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.qual199.txt

#-----------------------------------------  Prepare a combined Table  -------------------------------------------#
NFY_intersection_rep1_2_peaks <- read.table("NFY_intersection_rep1_2_peaks.bed", header = F, stringsAsFactors = F)
head(NFY_intersection_rep1_2_peaks)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6,NFY_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
tail(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- data.frame(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr[,c(21:23,1,24)])
colnames(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr) <- c("V1","V2","V3","V4","V5")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- data.frame(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr[order(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr$V4),]
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
head(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
tail(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr)

AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq <- data.frame(unique(AlspNFYddsNormcounts_B6_monofilt_uniqB6$V4))
colnames(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq) <- c("V4")
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq)

NFY_B6_specific_peaks_CCAAT <- read.table("NFY_B6_specific_peaks_CCAAT.bed")
NFY_B6_specific_peaks_CCAAT_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_CCAAT$V4),unique(NFY_B6_specific_peaks_CCAAT$V4))
colnames(NFY_B6_specific_peaks_CCAAT_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_CCAAT_uniq <- data.frame(NFY_B6_specific_peaks_CCAAT_uniq)
NFY_B6_specific_peaks_CCAAT_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_CCAAT_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_CCAAT_uniqB6 <- NFY_B6_specific_peaks_CCAAT_uniqB6[order(NFY_B6_specific_peaks_CCAAT_uniqB6$V4),]
dim(NFY_B6_specific_peaks_CCAAT_uniqB6)


NFY_B6_specific_peaks_jf1v2 <- read.table("NFY_B6_specific_peaks_jf1v2.bed")
NFY_B6_specific_peaks_jf1v2_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_jf1v2$V4),unique(NFY_B6_specific_peaks_jf1v2$V4))
colnames(NFY_B6_specific_peaks_jf1v2_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_jf1v2_uniq <- data.frame(NFY_B6_specific_peaks_jf1v2_uniq)
NFY_B6_specific_peaks_jf1v2_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_jf1v2_uniqB6 <- NFY_B6_specific_peaks_jf1v2_uniqB6[order(NFY_B6_specific_peaks_jf1v2_uniqB6$V4),]
dim(NFY_B6_specific_peaks_jf1v2_uniqB6)


NFY_B6_specific_peaks_CCAAT_jf1v2 <- read.table("NFY_B6_specific_peaks_CCAAT_jf1v2.bed")
NFY_B6_specific_peaks_CCAAT_jf1v2_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_CCAAT_jf1v2$V4),unique(NFY_B6_specific_peaks_CCAAT_jf1v2$V4))
colnames(NFY_B6_specific_peaks_CCAAT_jf1v2_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_CCAAT_jf1v2_uniq <- data.frame(NFY_B6_specific_peaks_CCAAT_jf1v2_uniq)
NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_CCAAT_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6 <- NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6[order(NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6$V4),]
dim(NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6)

NFY_B6_specific_peaks_CCAAT_then_jf1v2 <- read.table("NFY_B6_specific_peaks_CCAAT_then_jf1v2.bed")
NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_CCAAT_then_jf1v2$V4),unique(NFY_B6_specific_peaks_CCAAT_then_jf1v2$V4))
colnames(NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniq <- data.frame(NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniq)
NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6 <- NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6[order(NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6$V4),]
dim(NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6)#No need to include


Nomotif_NFY_B6_specific_peaks <- read.table("Nomotif_NFY_B6_specific_peaks.bed")
Nomotif_NFY_B6_specific_peaks_uniq <- cbind.data.frame(unique(Nomotif_NFY_B6_specific_peaks$V4),unique(Nomotif_NFY_B6_specific_peaks$V4))
colnames(Nomotif_NFY_B6_specific_peaks_uniq) <- c("V1","V4")
Nomotif_NFY_B6_specific_peaks_uniq <- data.frame(Nomotif_NFY_B6_specific_peaks_uniq)
Nomotif_NFY_B6_specific_peaks_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,Nomotif_NFY_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
Nomotif_NFY_B6_specific_peaks_uniqB6 <- Nomotif_NFY_B6_specific_peaks_uniqB6[order(Nomotif_NFY_B6_specific_peaks_uniqB6$V4),]
dim(Nomotif_NFY_B6_specific_peaks_uniqB6)


NomotifSNP_NFY_B6_specific_peaks <- read.table("NomotifSNP_NFY_B6_specific_peaks.bed")
NomotifSNP_NFY_B6_specific_peaks_uniq <- cbind.data.frame(unique(NomotifSNP_NFY_B6_specific_peaks$V4),unique(NomotifSNP_NFY_B6_specific_peaks$V4))
colnames(NomotifSNP_NFY_B6_specific_peaks_uniq) <- c("V1","V4")
NomotifSNP_NFY_B6_specific_peaks_uniq <- data.frame(NomotifSNP_NFY_B6_specific_peaks_uniq)
NomotifSNP_NFY_B6_specific_peaks_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NomotifSNP_NFY_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
NomotifSNP_NFY_B6_specific_peaks_uniqB6 <- NomotifSNP_NFY_B6_specific_peaks_uniqB6[order(NomotifSNP_NFY_B6_specific_peaks_uniqB6$V4),]
dim(NomotifSNP_NFY_B6_specific_peaks_uniqB6)

NoSNP_NFY_B6_specific_peaks <- read.table("NoSNP_NFY_B6_specific_peaks.bed")
NoSNP_NFY_B6_specific_peaks_uniq <- cbind.data.frame(unique(NoSNP_NFY_B6_specific_peaks$V4),unique(NoSNP_NFY_B6_specific_peaks$V4))
colnames(NoSNP_NFY_B6_specific_peaks_uniq) <- c("V1","V4")
NoSNP_NFY_B6_specific_peaks_uniq <- data.frame(NoSNP_NFY_B6_specific_peaks_uniq)
NoSNP_NFY_B6_specific_peaks_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NoSNP_NFY_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
NoSNP_NFY_B6_specific_peaks_uniqB6 <- NoSNP_NFY_B6_specific_peaks_uniqB6[order(NoSNP_NFY_B6_specific_peaks_uniqB6$V4),]
dim(NoSNP_NFY_B6_specific_peaks_uniqB6)


AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndel <- read.table("AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndel.txt")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu <- cbind.data.frame(unique(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndel$V4),unique(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndel$V4))
colnames(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu) <- c("V1","V4")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu <- data.frame(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq[order(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq$V4),]
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq)

NFY_B6_specific_peaks_CCAAT_jf1vindel <- read.table("NFY_B6_specific_peaks_CCAAT_jf1vindel.bed")
NFY_B6_specific_peaks_CCAAT_jf1vindel_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_CCAAT_jf1vindel$V4),unique(NFY_B6_specific_peaks_CCAAT_jf1vindel$V4))
colnames(NFY_B6_specific_peaks_CCAAT_jf1vindel_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_CCAAT_jf1vindel_uniq <- data.frame(NFY_B6_specific_peaks_CCAAT_jf1vindel_uniq)
NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_CCAAT_jf1vindel_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6 <- NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6[order(NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6$V4),]
dim(NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6)

AlspNFYddsNormcounts_B6_BJ1InputBi <- read.table("AlspNFYddsNormcounts_B6_BJ1InputBi.txt")
AlspNFYddsNormcounts_B6_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_B6_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_B6_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_B6_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_B6_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_B6_BJ1InputBiu)
AlspNFYddsNormcounts_B6_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_BJ1InputBiuniq <- AlspNFYddsNormcounts_B6_BJ1InputBiuniq[order(AlspNFYddsNormcounts_B6_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_B6_BJ1InputBiuniq)

AlspNFYddsNormcounts_B6_InputBiallelicfilt <- read.table("AlspNFYddsNormcounts_B6_InputBiallelicfilt.txt")
AlspNFYddsNormcounts_B6_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNFYddsNormcounts_B6_InputBiallelicfilt$V4),unique(AlspNFYddsNormcounts_B6_InputBiallelicfilt$V4))
colnames(AlspNFYddsNormcounts_B6_InputBiallelicfiltu) <- c("V1","V4")
AlspNFYddsNormcounts_B6_InputBiallelicfiltu <- data.frame(AlspNFYddsNormcounts_B6_InputBiallelicfiltu)
AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq <- AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq[order(AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq$V4),]
dim(AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq)

AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi <- read.table("AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi.txt")
AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiu)
AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq <- AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq[order(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq)

AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR <- read.table("AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu <- cbind.data.frame(unique(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR$V4),unique(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICR$V4))
colnames(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu) <- c("V1","V4")
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu <- data.frame(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq[order(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V4),]
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq)

#Its a blank file but needed, So 
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V1 <- NA
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq[order(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V4),]
dim(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq)

Alsp_NFY_B6_specific_peaks_CGIs_mm10_chr <- read.table("Alsp_NFY_B6_specific_peaks_CGIs_mm10_chr.txt")
Alsp_NFY_B6_specific_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chr$V4),unique(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NFY_B6_specific_peaks_CGIs_mm10_chru <- data.frame(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chru)
Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, Alsp_NFY_B6_specific_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq <- Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq[order(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq)

NFY_B6_specific_peaks_SP1_motif_SNP_Indels <- read.table("NFY_B6_specific_peaks_SP1_motif_SNP_Indels.bed")
NFY_B6_specific_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_B6_specific_peaks_SP1_motif_SNP_Indels$V4),unique(NFY_B6_specific_peaks_SP1_motif_SNP_Indels$V4))
colnames(NFY_B6_specific_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_B6_specific_peaks_SP1_motif_SNP_Indelsu <- data.frame(NFY_B6_specific_peaks_SP1_motif_SNP_Indelsu)
NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, NFY_B6_specific_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq[order(NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq)

NFY_B6_specific_peaks_YY1_motif_SNP_Indels <- read.table("NFY_B6_specific_peaks_YY1_motif_SNP_Indels.bed")
NFY_B6_specific_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_B6_specific_peaks_YY1_motif_SNP_Indels$V4),unique(NFY_B6_specific_peaks_YY1_motif_SNP_Indels$V4))
colnames(NFY_B6_specific_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_B6_specific_peaks_YY1_motif_SNP_Indelsu <- data.frame(NFY_B6_specific_peaks_YY1_motif_SNP_Indelsu)
NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, NFY_B6_specific_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq <- NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq[order(NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq)

NFY_B6_specific_peaks_NRF1_motif_SNP_Indels <- read.table("NFY_B6_specific_peaks_NRF1_motif_SNP_Indels.bed")
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_B6_specific_peaks_NRF1_motif_SNP_Indels$V4),unique(NFY_B6_specific_peaks_NRF1_motif_SNP_Indels$V4))
colnames(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsu <- data.frame(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsu)
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq)

#Its a blank file but needed, So 
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq <- AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq$V1 <- NA
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq)


#Pub NFY_B6
NFY_B6_specific_peaks_pub_NFY <- read.table("NFY_B6_specific_peaks_pub_NFY.bed")
NFY_B6_specific_peaks_pub_NFY_uniq <- cbind.data.frame(unique(NFY_B6_specific_peaks_pub_NFY$V4),unique(NFY_B6_specific_peaks_pub_NFY$V4))
colnames(NFY_B6_specific_peaks_pub_NFY_uniq) <- c("V1","V4")
NFY_B6_specific_peaks_pub_NFY_uniq <- data.frame(NFY_B6_specific_peaks_pub_NFY_uniq)
NFY_B6_specific_peaks_pub_NFY_uniqB6 <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq,NFY_B6_specific_peaks_pub_NFY_uniq, by="V4", all=T,sort = F)
NFY_B6_specific_peaks_pub_NFY_uniqB6 <- NFY_B6_specific_peaks_pub_NFY_uniqB6[order(NFY_B6_specific_peaks_pub_NFY_uniqB6$V4),]
dim(NFY_B6_specific_peaks_pub_NFY_uniqB6)

#Position Het SNPs
AlspNFYddsNormcounts_B6_BJ1InputBi <- read.table("AlspNFYddsNormcounts_B6_BJ1InputBi.txt")
AlspNFYddsNormcounts_B6_BJ1InputBiSNP <- AlspNFYddsNormcounts_B6_BJ1InputBi[,c(4,21:27)]
head(AlspNFYddsNormcounts_B6_BJ1InputBiSNP)
AlspNFYddsNormcounts_B6_BJ1InputBiSNP <- data.frame(AlspNFYddsNormcounts_B6_BJ1InputBiSNP)
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq <- merge(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNFYddsNormcounts_B6_BJ1InputBiSNP, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq["SNP"] <- paste0(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V21,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V22,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V23,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V24,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V25,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V26,"_",
                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq <- data.frame(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)

#Assign Score 0,1
NFY_B6_specific_peaks_CCAAT_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_CCAAT_uniqB6$V1),0,1)
NFY_B6_specific_peaks_jf1v2_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_jf1v2_uniqB6$V1),0,1)
NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6$V1),0,1)
NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_CCAAT_then_jf1v2_uniqB6$V1),0,1)
Nomotif_NFY_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(Nomotif_NFY_B6_specific_peaks_uniqB6$V1),0,1)
NomotifSNP_NFY_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(NomotifSNP_NFY_B6_specific_peaks_uniqB6$V1),0,1)
NoSNP_NFY_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(NoSNP_NFY_B6_specific_peaks_uniqB6$V1),0,1)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq$V1),0,1)
NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6$V1),0,1)
AlspNFYddsNormcounts_B6_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_B6_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq$V1),0,1)
AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V1),0,1)
Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq$V1),0,1)
NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq$V1),0,1)
NFY_B6_specific_peaks_pub_NFY_uniqB6["score"] <- ifelse(is.na(NFY_B6_specific_peaks_pub_NFY_uniqB6$V1),0,1)

matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features <- cbind.data.frame(AlspNFYddsNormcounts_B6_monofilt_uniqB6_uniq_chr,
                                                                            NFY_B6_specific_peaks_CCAAT_uniqB6,
                                                                            NFY_B6_specific_peaks_jf1v2_uniqB6,
                                                                            NFY_B6_specific_peaks_CCAAT_jf1v2_uniqB6,
                                                                            AlspNFYddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq,
                                                                            NFY_B6_specific_peaks_CCAAT_jf1vindel_uniqB6,
                                                                            AlspNFYddsNormcounts_B6_BJ1InputBiuniq,
                                                                            AlspNFYddsNormcounts_B6_InputBiallelicfiltuniq,
                                                                            AlspNFYddsNormcounts_B6_InputBiallelicfilt_BJ1InputBiuniq,
                                                                            AlspNFYddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq,
                                                                            Alsp_NFY_B6_specific_peaks_CGIs_mm10_chruniq,
                                                                            NFY_B6_specific_peaks_SP1_motif_SNP_Indelsuniq,
                                                                            NFY_B6_specific_peaks_YY1_motif_SNP_Indelsuniq,
                                                                            NFY_B6_specific_peaks_NRF1_motif_SNP_Indelsuniq,
                                                                            NFY_B6_specific_peaks_pub_NFY_uniqB6,
                                                                            AlspNFYddsNormcounts_B6_BJ1InputBiSNPuniq)

colnames(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features) <- c("Chr","Start","End","Peaks","Score","IDX.1","IDY.1","CCAAT","IDX.2","IDY.2","SNPs","IDX.3","IDY.3","CCAAT_SNPs","IDX.4","IDY.4","Indels","IDX.5","IDY.5","CCAAT_Indels","IDX.6","IDY.6","Input_Het_SNPs","IDX.7","IDY.7","Input_Bi_Regions","IDX.8","IDY.8","Input_Bi_Het_SNPs","IDX.9","IDY.9","ICRs","IDX.10","IDY.10","CGIs","IDX.11","IDY.11","SP1","IDX.12","IDY.12","YY1","IDX.13","IDY.13","NRF1","IDX.14","IDY.14","pub_NFY","IDX.15","PosHetSNPs")

dim(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features)
head(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features)
tail(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features)
matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab  <- data.frame(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features[,c(1:4,seq(5,47,3),49)])
matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab["Support_Index"] <- rowMeans(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab[,6:(length(colnames(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab))-1)])
matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab <- matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab[order(-matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab)
dim(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab)
write.table(matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab, "matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)


#-------- NFY_JF1
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1,NFY_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
dim(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
tail(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- data.frame(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr[,c(21:23,1,24)])
colnames(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr) <- c("V1","V2","V3","V4","V5")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- data.frame(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr[order(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr$V4),]
dim(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
head(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
tail(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)

AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq <- data.frame(unique(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1$V4))
colnames(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq) <- c("V4")
dim(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq)

NFY_JF1_specific_peaks_CCAAT <- read.table("NFY_JF1_specific_peaks_CCAAT.bed")
NFY_JF1_specific_peaks_CCAAT_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_CCAAT$V4),unique(NFY_JF1_specific_peaks_CCAAT$V4))
colnames(NFY_JF1_specific_peaks_CCAAT_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_CCAAT_uniq <- data.frame(NFY_JF1_specific_peaks_CCAAT_uniq)
NFY_JF1_specific_peaks_CCAAT_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_CCAAT_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_CCAAT_uniqJF1 <- NFY_JF1_specific_peaks_CCAAT_uniqJF1[order(NFY_JF1_specific_peaks_CCAAT_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_CCAAT_uniqJF1)


NFY_JF1_specific_peaks_jf1v2 <- read.table("NFY_JF1_specific_peaks_jf1v2.bed")
NFY_JF1_specific_peaks_jf1v2_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_jf1v2$V4),unique(NFY_JF1_specific_peaks_jf1v2$V4))
colnames(NFY_JF1_specific_peaks_jf1v2_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_jf1v2_uniq <- data.frame(NFY_JF1_specific_peaks_jf1v2_uniq)
NFY_JF1_specific_peaks_jf1v2_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_jf1v2_uniqJF1 <- NFY_JF1_specific_peaks_jf1v2_uniqJF1[order(NFY_JF1_specific_peaks_jf1v2_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_jf1v2_uniqJF1)


NFY_JF1_specific_peaks_CCAAT_jf1v2 <- read.table("NFY_JF1_specific_peaks_CCAAT_jf1v2.bed")
NFY_JF1_specific_peaks_CCAAT_jf1v2_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_CCAAT_jf1v2$V4),unique(NFY_JF1_specific_peaks_CCAAT_jf1v2$V4))
colnames(NFY_JF1_specific_peaks_CCAAT_jf1v2_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_CCAAT_jf1v2_uniq <- data.frame(NFY_JF1_specific_peaks_CCAAT_jf1v2_uniq)
NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_CCAAT_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1 <- NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1[order(NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1)

NFY_JF1_specific_peaks_CCAAT_then_jf1v2 <- read.table("NFY_JF1_specific_peaks_CCAAT_then_jf1v2.bed")
NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_CCAAT_then_jf1v2$V4),unique(NFY_JF1_specific_peaks_CCAAT_then_jf1v2$V4))
colnames(NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniq <- data.frame(NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniq)
NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1 <- NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1[order(NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1)#No need to include


Nomotif_NFY_JF1_specific_peaks <- read.table("Nomotif_NFY_JF1_specific_peaks.bed")
Nomotif_NFY_JF1_specific_peaks_uniq <- cbind.data.frame(unique(Nomotif_NFY_JF1_specific_peaks$V4),unique(Nomotif_NFY_JF1_specific_peaks$V4))
colnames(Nomotif_NFY_JF1_specific_peaks_uniq) <- c("V1","V4")
Nomotif_NFY_JF1_specific_peaks_uniq <- data.frame(Nomotif_NFY_JF1_specific_peaks_uniq)
Nomotif_NFY_JF1_specific_peaks_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,Nomotif_NFY_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
Nomotif_NFY_JF1_specific_peaks_uniqJF1 <- Nomotif_NFY_JF1_specific_peaks_uniqJF1[order(Nomotif_NFY_JF1_specific_peaks_uniqJF1$V4),]
dim(Nomotif_NFY_JF1_specific_peaks_uniqJF1)


NomotifSNP_NFY_JF1_specific_peaks <- read.table("NomotifSNP_NFY_JF1_specific_peaks.bed")
NomotifSNP_NFY_JF1_specific_peaks_uniq <- cbind.data.frame(unique(NomotifSNP_NFY_JF1_specific_peaks$V4),unique(NomotifSNP_NFY_JF1_specific_peaks$V4))
colnames(NomotifSNP_NFY_JF1_specific_peaks_uniq) <- c("V1","V4")
NomotifSNP_NFY_JF1_specific_peaks_uniq <- data.frame(NomotifSNP_NFY_JF1_specific_peaks_uniq)
NomotifSNP_NFY_JF1_specific_peaks_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NomotifSNP_NFY_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
NomotifSNP_NFY_JF1_specific_peaks_uniqJF1 <- NomotifSNP_NFY_JF1_specific_peaks_uniqJF1[order(NomotifSNP_NFY_JF1_specific_peaks_uniqJF1$V4),]
dim(NomotifSNP_NFY_JF1_specific_peaks_uniqJF1)

NoSNP_NFY_JF1_specific_peaks <- read.table("NoSNP_NFY_JF1_specific_peaks.bed")
NoSNP_NFY_JF1_specific_peaks_uniq <- cbind.data.frame(unique(NoSNP_NFY_JF1_specific_peaks$V4),unique(NoSNP_NFY_JF1_specific_peaks$V4))
colnames(NoSNP_NFY_JF1_specific_peaks_uniq) <- c("V1","V4")
NoSNP_NFY_JF1_specific_peaks_uniq <- data.frame(NoSNP_NFY_JF1_specific_peaks_uniq)
NoSNP_NFY_JF1_specific_peaks_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NoSNP_NFY_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
NoSNP_NFY_JF1_specific_peaks_uniqJF1 <- NoSNP_NFY_JF1_specific_peaks_uniqJF1[order(NoSNP_NFY_JF1_specific_peaks_uniqJF1$V4),]
dim(NoSNP_NFY_JF1_specific_peaks_uniqJF1)

AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel <- read.table("AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel.txt")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu <- cbind.data.frame(unique(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel$V4),unique(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel$V4))
colnames(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu) <- c("V1","V4")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu <- data.frame(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq[order(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq$V4),]
dim(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq)

NFY_JF1_specific_peaks_CCAAT_jf1vindel <- read.table("NFY_JF1_specific_peaks_CCAAT_jf1vindel.bed")
NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_CCAAT_jf1vindel$V4),unique(NFY_JF1_specific_peaks_CCAAT_jf1vindel$V4))
colnames(NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniq <- data.frame(NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniq)
NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1 <- NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1[order(NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1)

AlspNFYddsNormcounts_JF1_BJ1InputBi <- read.table("AlspNFYddsNormcounts_JF1_BJ1InputBi.txt")
AlspNFYddsNormcounts_JF1_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_JF1_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_JF1_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_JF1_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_JF1_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_JF1_BJ1InputBiu)
AlspNFYddsNormcounts_JF1_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_BJ1InputBiuniq <- AlspNFYddsNormcounts_JF1_BJ1InputBiuniq[order(AlspNFYddsNormcounts_JF1_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_BJ1InputBiuniq)

AlspNFYddsNormcounts_JF1_InputBiallelicfilt <- read.table("AlspNFYddsNormcounts_JF1_InputBiallelicfilt.txt")
AlspNFYddsNormcounts_JF1_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNFYddsNormcounts_JF1_InputBiallelicfilt$V4),unique(AlspNFYddsNormcounts_JF1_InputBiallelicfilt$V4))
colnames(AlspNFYddsNormcounts_JF1_InputBiallelicfiltu) <- c("V1","V4")
AlspNFYddsNormcounts_JF1_InputBiallelicfiltu <- data.frame(AlspNFYddsNormcounts_JF1_InputBiallelicfiltu)
AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq <- AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq[order(AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq)

AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi <- read.table("AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi.txt")
AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiu)
AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq <- AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq[order(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq)

AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR <- read.table("AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu <- cbind.data.frame(unique(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR$V4),unique(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR$V4))
colnames(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu) <- c("V1","V4")
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu <- data.frame(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq[order(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq)


Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chr <- read.table("Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chr.txt")
Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chr$V4),unique(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chru <- data.frame(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chru)
Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq <- Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq[order(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq)

NFY_JF1_specific_peaks_SP1_motif_SNP_Indels <- read.table("NFY_JF1_specific_peaks_SP1_motif_SNP_Indels.bed")
NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_JF1_specific_peaks_SP1_motif_SNP_Indels$V4),unique(NFY_JF1_specific_peaks_SP1_motif_SNP_Indels$V4))
colnames(NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsu <- data.frame(NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsu)
NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq <- NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq[order(NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq)

NFY_JF1_specific_peaks_YY1_motif_SNP_Indels <- read.table("NFY_JF1_specific_peaks_YY1_motif_SNP_Indels.bed")
NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_JF1_specific_peaks_YY1_motif_SNP_Indels$V4),unique(NFY_JF1_specific_peaks_YY1_motif_SNP_Indels$V4))
colnames(NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsu <- data.frame(NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsu)
NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq <- NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq[order(NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq)

NFY_JF1_specific_peaks_NRF1_motif_SNP_Indels <- read.table("NFY_JF1_specific_peaks_NRF1_motif_SNP_Indels.bed")
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indels$V4),unique(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indels$V4))
colnames(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsu <- data.frame(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsu)
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq)

#Its a blank file but needed, So 
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq <- AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq$V1 <- NA
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq)

#Pub NFY_JF1
NFY_JF1_specific_peaks_pub_NFY <- read.table("NFY_JF1_specific_peaks_pub_NFY.bed")
NFY_JF1_specific_peaks_pub_NFY_uniq <- cbind.data.frame(unique(NFY_JF1_specific_peaks_pub_NFY$V4),unique(NFY_JF1_specific_peaks_pub_NFY$V4))
colnames(NFY_JF1_specific_peaks_pub_NFY_uniq) <- c("V1","V4")
NFY_JF1_specific_peaks_pub_NFY_uniq <- data.frame(NFY_JF1_specific_peaks_pub_NFY_uniq)
NFY_JF1_specific_peaks_pub_NFY_uniqJF1 <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq,NFY_JF1_specific_peaks_pub_NFY_uniq, by="V4", all=T,sort = F)
NFY_JF1_specific_peaks_pub_NFY_uniqJF1 <- NFY_JF1_specific_peaks_pub_NFY_uniqJF1[order(NFY_JF1_specific_peaks_pub_NFY_uniqJF1$V4),]
dim(NFY_JF1_specific_peaks_pub_NFY_uniqJF1)

#Position Het SNPs
AlspNFYddsNormcounts_JF1_BJ1InputBi <- read.table("AlspNFYddsNormcounts_JF1_BJ1InputBi.txt")
AlspNFYddsNormcounts_JF1_BJ1InputBiSNP <- AlspNFYddsNormcounts_JF1_BJ1InputBi[,c(4,21:27)]
head(AlspNFYddsNormcounts_JF1_BJ1InputBiSNP)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNP <- data.frame(AlspNFYddsNormcounts_JF1_BJ1InputBiSNP)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq <- merge(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNFYddsNormcounts_JF1_BJ1InputBiSNP, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq["SNP"] <- paste0(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V21,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V22,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V23,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V24,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V25,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V26,"_",
                                                             AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq <- data.frame(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)

NFY_JF1_specific_peaks_CCAAT_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_CCAAT_uniqJF1$V1),0,1)
NFY_JF1_specific_peaks_jf1v2_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_jf1v2_uniqJF1$V1),0,1)
NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1$V1),0,1)
NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_CCAAT_then_jf1v2_uniqJF1$V1),0,1)
Nomotif_NFY_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(Nomotif_NFY_JF1_specific_peaks_uniqJF1$V1),0,1)
NomotifSNP_NFY_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(NomotifSNP_NFY_JF1_specific_peaks_uniqJF1$V1),0,1)
NoSNP_NFY_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(NoSNP_NFY_JF1_specific_peaks_uniqJF1$V1),0,1)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq$V1),0,1)
NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1$V1),0,1)
AlspNFYddsNormcounts_JF1_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_JF1_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq$V1),0,1)
AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V1),0,1)
Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq$V1),0,1)
NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq$V1),0,1)
NFY_JF1_specific_peaks_pub_NFY_uniqJF1["score"] <- ifelse(is.na(NFY_JF1_specific_peaks_pub_NFY_uniqJF1$V1),0,1)


matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features <- cbind.data.frame(AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr,
                                                                              NFY_JF1_specific_peaks_CCAAT_uniqJF1,
                                                                              NFY_JF1_specific_peaks_jf1v2_uniqJF1,
                                                                              NFY_JF1_specific_peaks_CCAAT_jf1v2_uniqJF1,
                                                                              AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq,
                                                                              NFY_JF1_specific_peaks_CCAAT_jf1vindel_uniqJF1,
                                                                              AlspNFYddsNormcounts_JF1_BJ1InputBiuniq,
                                                                              AlspNFYddsNormcounts_JF1_InputBiallelicfiltuniq,
                                                                              AlspNFYddsNormcounts_JF1_InputBiallelicfilt_BJ1InputBiuniq,
                                                                              AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq,
                                                                              Alsp_NFY_JF1_specific_peaks_CGIs_mm10_chruniq,
                                                                              NFY_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq,
                                                                              NFY_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq,
                                                                              NFY_JF1_specific_peaks_NRF1_motif_SNP_Indelsuniq,
                                                                              NFY_JF1_specific_peaks_pub_NFY_uniqJF1,
                                                                              AlspNFYddsNormcounts_JF1_BJ1InputBiSNPuniq)



colnames(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features) <- c("Chr","Start","End","Peaks","Score","IDX.1","IDY.1","CCAAT","IDX.2","IDY.2","SNPs","IDX.3","IDY.3","CCAAT_SNPs","IDX.4","IDY.4","Indels","IDX.5","IDY.5","CCAAT_Indels","IDX.6","IDY.6","Input_Het_SNPs","IDX.7","IDY.7","Input_Bi_Regions","IDX.8","IDY.8","Input_Bi_Het_SNPs","IDX.9","IDY.9","ICRs","IDX.10","IDY.10","CGIs","IDX.11","IDY.11","SP1","IDX.12","IDY.12","YY1","IDX.13","IDY.13","NRF1","IDX.14","IDY.14","pub_NFY","IDX.15","PosHetSNPs")

dim(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features)
head(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features)
tail(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features)
matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab  <- data.frame(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features[,c(1:4,seq(5,47,3),49)])
matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab["Support_Index"] <- rowMeans(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab[,6:(length(colnames(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab))-1)])
matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab <- matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab[order(-matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab)

write.table(matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab, "matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)


#----------------------------- Biallelic peaks intersection ------------------
#AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt

#Biallelic
#Biallelicalleleic
AlspInputddsNormcounts_Biallelic <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio <= 0.60 &
                                                                                     AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio >= 0.40),]
head(AlspInputddsNormcounts_Biallelic)
dim(AlspInputddsNormcounts_Biallelic)

#Bialleically expressed  with FDR > 0.05
AlspInputddsNormcounts_Biallelicfilt <- AlspInputddsNormcounts_Biallelic[which(AlspInputddsNormcounts_Biallelic$Input_alavg.comb.fdr>0.05),]
head(AlspInputddsNormcounts_Biallelicfilt)
dim(AlspInputddsNormcounts_Biallelicfilt)
write.table(AlspInputddsNormcounts_Biallelicfilt,"AlspInputddsNormcounts_Biallelicfilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)



#NFY
#Motif and SNP overlap peaks
#CCAAT
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_motif.mm10.sort.bed > NFY_B_biallelic_peaks_CCAAT.bed
#SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_B_biallelic_peaks_jf1v2.bed

#Note CCAAT_jf1v2SNPmm10.bed is file intersected with SNP (jf1v2mm10)
#CCAAT+SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1v2SNPmm10.bed > NFY_B_biallelic_peaks_CCAAT_jf1v2.bed

#CCAAT, SNP
bedtools intersect -wa -wb -a NFY_B_biallelic_peaks_CCAAT.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_B_biallelic_peaks_CCAAT_then_jf1v2.bed

#NO Nfy
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_motif.mm10.sort.bed -v > Nomotif_NFY_B_biallelic_peaks.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1v2SNPmm10.bed -v > NomotifSNP_NFY_B_biallelic_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NFY_B_biallelic_peaks.bed

#Motif and SNP overlap peaks JF1
#CCAAT
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_motif.mm10.jf1.sort.bed > NFY_J_biallelic_peaks_CCAAT.bed
#SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_J_biallelic_peaks_jf1v2.bed
#Note CCAAT_jf1v2SNPmm10.jf1.bed is file intersected with SNP (jf1v2mm10.jf1)
#CCAAT+SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1v2SNPmm10.jf1.bed > NFY_J_biallelic_peaks_CCAAT_jf1v2.bed
#CCAAT, SNP
bedtools intersect -wa -wb -a NFY_J_biallelic_peaks_CCAAT.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NFY_J_biallelic_peaks_CCAAT_then_jf1v2.bed

#NO Nfy
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_motif.mm10.jf1.sort.bed -v > Nomotif_NFY_J_biallelic_peaks.bed
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1v2SNPmm10.jf1.bed -v > NomotifSNP_NFY_J_biallelic_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NFY_J_biallelic_peaks.bed


bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b BJ1resultInputBi_chr.txt > AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b AlspInputddsNormcounts_Biallelicfilt_BJ1resultInputBi.txt > AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi.txt
bedtools intersect -wa -wb -a  AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt  -b AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt > AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt

#Indel
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt

#Motif, Indels, as before
#CCAAT+Indels,B6 based
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1vindelSNPmm10.bed > NFY_B_biallelic_peaks_CCAAT_jf1vindel.bed

#CCAAT+Indels,JF1-based
bedtools intersect -wa -wb -a AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b CCAAT_jf1vindelSNPmm10.jf1.bed > NFY_J_biallelic_peaks_CCAAT_jf1vindel.bed


#Countings
#Common NFY Biallelic B and J
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_qual199.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt | wc -l
sort -k4,4 -u AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt | wc -l

#NFY_B
sort -k4,4 -u NFY_B_biallelic_peaks_CCAAT_jf1vindel.bed | wc -l
sort -k4,4 -u NFY_B_biallelic_peaks_CCAAT.bed | wc -l
sort -k4,4 -u NFY_B_biallelic_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NFY_B_biallelic_peaks_CCAAT_jf1v2.bed | wc -l
sort -k4,4 -u NFY_B_biallelic_peaks_CCAAT_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NFY_B_biallelic_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NFY_B_biallelic_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NFY_B_biallelic_peaks.bed | wc -l

#NFY_J
sort -k4,4 -u NFY_J_biallelic_peaks_CCAAT_jf1vindel.bed | wc -l
sort -k4,4 -u NFY_J_biallelic_peaks_CCAAT.bed | wc -l
sort -k4,4 -u NFY_J_biallelic_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NFY_J_biallelic_peaks_CCAAT_jf1v2.bed | wc -l
sort -k4,4 -u NFY_J_biallelic_peaks_CCAAT_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NFY_J_biallelic_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NFY_J_biallelic_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NFY_J_biallelic_peaks.bed | wc -l


awk '{if($5>=199) print $0}' AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.txt > AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1.qual199.txt

##-- Prepare Table NFY Biallelic Peaks --##
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- merge(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1,NFY_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
dim(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
tail(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr[,c(21:23,1,24)])
colnames(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr) <- c("V1","V2","V3","V4","V5")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr[order(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
head(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
tail(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)

AlspNFYddsNormcounts_Biallelicfilt_shared <- data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_shared) <- c("V4")
dim(AlspNFYddsNormcounts_Biallelicfilt_shared)


NFY_B_biallelic_peaks_CCAAT <- read.table("NFY_B_biallelic_peaks_CCAAT.bed")
NFY_B_biallelic_peaks_CCAAT_shared <- cbind.data.frame(unique(NFY_B_biallelic_peaks_CCAAT$V4),unique(NFY_B_biallelic_peaks_CCAAT$V4))
colnames(NFY_B_biallelic_peaks_CCAAT_shared) <- c("V1","V4")
NFY_B_biallelic_peaks_CCAAT_shared <- data.frame(NFY_B_biallelic_peaks_CCAAT_shared)
NFY_B_biallelic_peaks_CCAAT_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_B_biallelic_peaks_CCAAT_shared, by="V4", all=T,sort = F)
NFY_B_biallelic_peaks_CCAAT_shared <- NFY_B_biallelic_peaks_CCAAT_shared[order(NFY_B_biallelic_peaks_CCAAT_shared$V4),]
dim(NFY_B_biallelic_peaks_CCAAT_shared)


NFY_B_biallelic_peaks_jf1v2 <- read.table("NFY_B_biallelic_peaks_jf1v2.bed")
NFY_B_biallelic_peaks_jf1v2_shared <- cbind.data.frame(unique(NFY_B_biallelic_peaks_jf1v2$V4),unique(NFY_B_biallelic_peaks_jf1v2$V4))
colnames(NFY_B_biallelic_peaks_jf1v2_shared) <- c("V1","V4")
NFY_B_biallelic_peaks_jf1v2_shared <- data.frame(NFY_B_biallelic_peaks_jf1v2_shared)
NFY_B_biallelic_peaks_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_B_biallelic_peaks_jf1v2_shared, by="V4", all=T,sort = F)
NFY_B_biallelic_peaks_jf1v2_shared <- NFY_B_biallelic_peaks_jf1v2_shared[order(NFY_B_biallelic_peaks_jf1v2_shared$V4),]
dim(NFY_B_biallelic_peaks_jf1v2_shared)


NFY_B_biallelic_peaks_CCAAT_jf1v2 <- read.table("NFY_B_biallelic_peaks_CCAAT_jf1v2.bed")
NFY_B_biallelic_peaks_CCAAT_jf1v2_shared <- cbind.data.frame(unique(NFY_B_biallelic_peaks_CCAAT_jf1v2$V4),unique(NFY_B_biallelic_peaks_CCAAT_jf1v2$V4))
colnames(NFY_B_biallelic_peaks_CCAAT_jf1v2_shared) <- c("V1","V4")
NFY_B_biallelic_peaks_CCAAT_jf1v2_shared <- data.frame(NFY_B_biallelic_peaks_CCAAT_jf1v2_shared)
NFY_B_biallelic_peaks_CCAAT_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_B_biallelic_peaks_CCAAT_jf1v2_shared, by="V4", all=T,sort = F)
NFY_B_biallelic_peaks_CCAAT_jf1v2_shared <- NFY_B_biallelic_peaks_CCAAT_jf1v2_shared[order(NFY_B_biallelic_peaks_CCAAT_jf1v2_shared$V4),]
dim(NFY_B_biallelic_peaks_CCAAT_jf1v2_shared)

NFY_B_biallelic_peaks_CCAAT_then_jf1v2 <- read.table("NFY_B_biallelic_peaks_CCAAT_then_jf1v2.bed")
NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared <- cbind.data.frame(unique(NFY_B_biallelic_peaks_CCAAT_then_jf1v2$V4),unique(NFY_B_biallelic_peaks_CCAAT_then_jf1v2$V4))
colnames(NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared) <- c("V1","V4")
NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared <- data.frame(NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared)
NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared, by="V4", all=T,sort = F)
NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared <- NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared[order(NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared$V4),]
dim(NFY_B_biallelic_peaks_CCAAT_then_jf1v2_shared)#No need to include


Nomotif_NFY_B_biallelic_peaks <- read.table("Nomotif_NFY_B_biallelic_peaks.bed")
Nomotif_NFY_B_biallelic_peaks_shared <- cbind.data.frame(unique(Nomotif_NFY_B_biallelic_peaks$V4),unique(Nomotif_NFY_B_biallelic_peaks$V4))
colnames(Nomotif_NFY_B_biallelic_peaks_shared) <- c("V1","V4")
Nomotif_NFY_B_biallelic_peaks_shared <- data.frame(Nomotif_NFY_B_biallelic_peaks_shared)
Nomotif_NFY_B_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,Nomotif_NFY_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
Nomotif_NFY_B_biallelic_peaks_shared <- Nomotif_NFY_B_biallelic_peaks_shared[order(Nomotif_NFY_B_biallelic_peaks_shared$V4),]
dim(Nomotif_NFY_B_biallelic_peaks_shared)


NomotifSNP_NFY_B_biallelic_peaks <- read.table("NomotifSNP_NFY_B_biallelic_peaks.bed")
NomotifSNP_NFY_B_biallelic_peaks_shared <- cbind.data.frame(unique(NomotifSNP_NFY_B_biallelic_peaks$V4),unique(NomotifSNP_NFY_B_biallelic_peaks$V4))
colnames(NomotifSNP_NFY_B_biallelic_peaks_shared) <- c("V1","V4")
NomotifSNP_NFY_B_biallelic_peaks_shared <- data.frame(NomotifSNP_NFY_B_biallelic_peaks_shared)
NomotifSNP_NFY_B_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NomotifSNP_NFY_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
NomotifSNP_NFY_B_biallelic_peaks_shared <- NomotifSNP_NFY_B_biallelic_peaks_shared[order(NomotifSNP_NFY_B_biallelic_peaks_shared$V4),]
dim(NomotifSNP_NFY_B_biallelic_peaks_shared)

NoSNP_NFY_B_biallelic_peaks <- read.table("NoSNP_NFY_B_biallelic_peaks.bed")
NoSNP_NFY_B_biallelic_peaks_shared <- cbind.data.frame(unique(NoSNP_NFY_B_biallelic_peaks$V4),unique(NoSNP_NFY_B_biallelic_peaks$V4))
colnames(NoSNP_NFY_B_biallelic_peaks_shared) <- c("V1","V4")
NoSNP_NFY_B_biallelic_peaks_shared <- data.frame(NoSNP_NFY_B_biallelic_peaks_shared)
NoSNP_NFY_B_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NoSNP_NFY_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
NoSNP_NFY_B_biallelic_peaks_shared <- NoSNP_NFY_B_biallelic_peaks_shared[order(NoSNP_NFY_B_biallelic_peaks_shared$V4),]
dim(NoSNP_NFY_B_biallelic_peaks_shared)

NFY_B_biallelic_peaks_CCAAT_jf1vindel <- read.table("NFY_B_biallelic_peaks_CCAAT_jf1vindel.bed")
NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniq <- cbind.data.frame(unique(NFY_B_biallelic_peaks_CCAAT_jf1vindel$V4),unique(NFY_B_biallelic_peaks_CCAAT_jf1vindel$V4))
colnames(NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniq) <- c("V1","V4")
NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniq <- data.frame(NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniq)
NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6 <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared,NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniq, by="V4", all=T,sort = F)
NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6 <- NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6[order(NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6$V4),]
dim(NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6)



NFY_J_biallelic_peaks_CCAAT <- read.table("NFY_J_biallelic_peaks_CCAAT.bed")
NFY_J_biallelic_peaks_CCAAT_shared <- cbind.data.frame(unique(NFY_J_biallelic_peaks_CCAAT$V4),unique(NFY_J_biallelic_peaks_CCAAT$V4))
colnames(NFY_J_biallelic_peaks_CCAAT_shared) <- c("V1","V4")
NFY_J_biallelic_peaks_CCAAT_shared <- data.frame(NFY_J_biallelic_peaks_CCAAT_shared)
NFY_J_biallelic_peaks_CCAAT_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_J_biallelic_peaks_CCAAT_shared, by="V4", all=T,sort = F)
NFY_J_biallelic_peaks_CCAAT_shared <- NFY_J_biallelic_peaks_CCAAT_shared[order(NFY_J_biallelic_peaks_CCAAT_shared$V4),]
dim(NFY_J_biallelic_peaks_CCAAT_shared)


NFY_J_biallelic_peaks_jf1v2 <- read.table("NFY_J_biallelic_peaks_jf1v2.bed")
NFY_J_biallelic_peaks_jf1v2_shared <- cbind.data.frame(unique(NFY_J_biallelic_peaks_jf1v2$V4),unique(NFY_J_biallelic_peaks_jf1v2$V4))
colnames(NFY_J_biallelic_peaks_jf1v2_shared) <- c("V1","V4")
NFY_J_biallelic_peaks_jf1v2_shared <- data.frame(NFY_J_biallelic_peaks_jf1v2_shared)
NFY_J_biallelic_peaks_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_J_biallelic_peaks_jf1v2_shared, by="V4", all=T,sort = F)
NFY_J_biallelic_peaks_jf1v2_shared <- NFY_J_biallelic_peaks_jf1v2_shared[order(NFY_J_biallelic_peaks_jf1v2_shared$V4),]
dim(NFY_J_biallelic_peaks_jf1v2_shared)


NFY_J_biallelic_peaks_CCAAT_jf1v2 <- read.table("NFY_J_biallelic_peaks_CCAAT_jf1v2.bed")
NFY_J_biallelic_peaks_CCAAT_jf1v2_shared <- cbind.data.frame(unique(NFY_J_biallelic_peaks_CCAAT_jf1v2$V4),unique(NFY_J_biallelic_peaks_CCAAT_jf1v2$V4))
colnames(NFY_J_biallelic_peaks_CCAAT_jf1v2_shared) <- c("V1","V4")
NFY_J_biallelic_peaks_CCAAT_jf1v2_shared <- data.frame(NFY_J_biallelic_peaks_CCAAT_jf1v2_shared)
NFY_J_biallelic_peaks_CCAAT_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_J_biallelic_peaks_CCAAT_jf1v2_shared, by="V4", all=T,sort = F)
NFY_J_biallelic_peaks_CCAAT_jf1v2_shared <- NFY_J_biallelic_peaks_CCAAT_jf1v2_shared[order(NFY_J_biallelic_peaks_CCAAT_jf1v2_shared$V4),]
dim(NFY_J_biallelic_peaks_CCAAT_jf1v2_shared)

NFY_J_biallelic_peaks_CCAAT_then_jf1v2 <- read.table("NFY_J_biallelic_peaks_CCAAT_then_jf1v2.bed")
NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared <- cbind.data.frame(unique(NFY_J_biallelic_peaks_CCAAT_then_jf1v2$V4),unique(NFY_J_biallelic_peaks_CCAAT_then_jf1v2$V4))
colnames(NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared) <- c("V1","V4")
NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared <- data.frame(NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared)
NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared, by="V4", all=T,sort = F)
NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared <- NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared[order(NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared$V4),]
dim(NFY_J_biallelic_peaks_CCAAT_then_jf1v2_shared)#No need to include


Nomotif_NFY_J_biallelic_peaks <- read.table("Nomotif_NFY_J_biallelic_peaks.bed")
Nomotif_NFY_J_biallelic_peaks_shared <- cbind.data.frame(unique(Nomotif_NFY_J_biallelic_peaks$V4),unique(Nomotif_NFY_J_biallelic_peaks$V4))
colnames(Nomotif_NFY_J_biallelic_peaks_shared) <- c("V1","V4")
Nomotif_NFY_J_biallelic_peaks_shared <- data.frame(Nomotif_NFY_J_biallelic_peaks_shared)
Nomotif_NFY_J_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,Nomotif_NFY_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
Nomotif_NFY_J_biallelic_peaks_shared <- Nomotif_NFY_J_biallelic_peaks_shared[order(Nomotif_NFY_J_biallelic_peaks_shared$V4),]
dim(Nomotif_NFY_J_biallelic_peaks_shared)


NomotifSNP_NFY_J_biallelic_peaks <- read.table("NomotifSNP_NFY_J_biallelic_peaks.bed")
NomotifSNP_NFY_J_biallelic_peaks_shared <- cbind.data.frame(unique(NomotifSNP_NFY_J_biallelic_peaks$V4),unique(NomotifSNP_NFY_J_biallelic_peaks$V4))
colnames(NomotifSNP_NFY_J_biallelic_peaks_shared) <- c("V1","V4")
NomotifSNP_NFY_J_biallelic_peaks_shared <- data.frame(NomotifSNP_NFY_J_biallelic_peaks_shared)
NomotifSNP_NFY_J_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NomotifSNP_NFY_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
NomotifSNP_NFY_J_biallelic_peaks_shared <- NomotifSNP_NFY_J_biallelic_peaks_shared[order(NomotifSNP_NFY_J_biallelic_peaks_shared$V4),]
dim(NomotifSNP_NFY_J_biallelic_peaks_shared)

NoSNP_NFY_J_biallelic_peaks <- read.table("NoSNP_NFY_J_biallelic_peaks.bed")
NoSNP_NFY_J_biallelic_peaks_shared <- cbind.data.frame(unique(NoSNP_NFY_J_biallelic_peaks$V4),unique(NoSNP_NFY_J_biallelic_peaks$V4))
colnames(NoSNP_NFY_J_biallelic_peaks_shared) <- c("V1","V4")
NoSNP_NFY_J_biallelic_peaks_shared <- data.frame(NoSNP_NFY_J_biallelic_peaks_shared)
NoSNP_NFY_J_biallelic_peaks_shared <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared ,NoSNP_NFY_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
NoSNP_NFY_J_biallelic_peaks_shared <- NoSNP_NFY_J_biallelic_peaks_shared[order(NoSNP_NFY_J_biallelic_peaks_shared$V4),]
dim(NoSNP_NFY_J_biallelic_peaks_shared)

AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel <- read.table("AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu <- cbind.data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel$V4),unique(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu) <- c("V1","V4")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq[order(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq)

NFY_J_biallelic_peaks_CCAAT_jf1vindel <- read.table("NFY_J_biallelic_peaks_CCAAT_jf1vindel.bed")
NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniq <- cbind.data.frame(unique(NFY_J_biallelic_peaks_CCAAT_jf1vindel$V4),unique(NFY_J_biallelic_peaks_CCAAT_jf1vindel$V4))
colnames(NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniq) <- c("V1","V4")
NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniq <- data.frame(NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniq)
NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1 <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared,NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniq, by="V4", all=T,sort = F)
NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1 <- NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1[order(NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1$V4),]
dim(NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1)

AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi <- read.table("AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi.txt")
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiu)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq <- AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq[order(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq)

AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt <- read.table("AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt")
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt$V4),unique(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltu) <- c("V1","V4")
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltu <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltu)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq <- AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq[order(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq)

AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi <- read.table("AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi.txt")
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiu <- cbind.data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi$V4),unique(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBi$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiu) <- c("V1","V4")
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiu <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiu)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq <- AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq[order(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq)

AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR <- read.table("AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu <- cbind.data.frame(unique(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR$V4),unique(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR$V4))
colnames(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu) <- c("V1","V4")
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq[order(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq)

#Its a blank file but needed, So 
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNFYddsNormcounts_Biallelicfilt_shared
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V1 <- NA
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq[order(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq)


Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr <- read.table("Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr.txt")
Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr$V4),unique(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru <- data.frame(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru)
Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq <- Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq[order(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq)


NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels <- read.table("NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels.bed")
NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels$V4),unique(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels$V4))
colnames(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu <- data.frame(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu)
NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq <- NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq[order(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq)

NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels <- read.table("NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels.bed")
NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels$V4),unique(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels$V4))
colnames(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu <- data.frame(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu)
NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq <- NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq[order(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq)


NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indels <- read.table("NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indels.bed")
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsu <- cbind.data.frame(unique(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indels$V4),unique(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indels$V4))
colnames(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsu) <- c("V1","V4")
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsu <- data.frame(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsu)
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq)

#Its a blank file but needed, So 
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq <- AlspNFYddsNormcounts_Biallelicfilt_shared
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq$V1 <- NA
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq <- NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq[order(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq$V4),]
dim(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq)

#Pub NFY_Biallelic
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY <- read.table("NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY.bed")
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniq <- cbind.data.frame(unique(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY$V4),unique(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY$V4))
colnames(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniq) <- c("V1","V4")
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniq <- data.frame(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniq)
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared,NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniq, by="V4", all=T,sort = F)
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi <- NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi[order(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi$V4),]
dim(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi)

#Position Het SNPs
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi <- read.table("AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi.txt")
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNP <- AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBi[,c(4,21:27)]
head(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNP)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNP <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNP)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq <- merge(AlspNFYddsNormcounts_Biallelicfilt_shared, AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNP, by="V4", all=T,sort = F)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq["SNP"] <- paste0(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V21,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V22,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V23,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V24,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V25,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V26,"_",
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq <- data.frame(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq <- AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq[order(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq$V4),]
dim(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)
head(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)

#Assign Score 0,1
NFY_B_biallelic_peaks_CCAAT_shared["score"] <- ifelse(is.na(NFY_B_biallelic_peaks_CCAAT_shared$V1),0,1)
NFY_B_biallelic_peaks_jf1v2_shared["score"] <- ifelse(is.na(NFY_B_biallelic_peaks_jf1v2_shared$V1),0,1)
NFY_B_biallelic_peaks_CCAAT_jf1v2_shared["score"] <- ifelse(is.na(NFY_B_biallelic_peaks_CCAAT_jf1v2_shared$V1),0,1)
NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6["score"] <- ifelse(is.na(NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6$V1),0,1)
NFY_J_biallelic_peaks_CCAAT_shared["score"] <- ifelse(is.na(NFY_J_biallelic_peaks_CCAAT_shared$V1),0,1)
NFY_J_biallelic_peaks_jf1v2_shared["score"] <- ifelse(is.na(NFY_J_biallelic_peaks_jf1v2_shared$V1),0,1)
NFY_J_biallelic_peaks_CCAAT_jf1v2_shared["score"] <- ifelse(is.na(NFY_J_biallelic_peaks_CCAAT_jf1v2_shared$V1),0,1)
NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1["score"] <- ifelse(is.na(NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1$V1),0,1)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq$V1),0,1)
AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq$V1),0,1)
AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq$V1),0,1)
AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V1),0,1)
Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq$V1),0,1)
NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq$V1),0,1)
NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi["score"] <- ifelse(is.na(NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi$V1),0,1)

matrix_AlspNFYddsNormcounts_Biallelicfilt_features <- cbind.data.frame(AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr,
                                                                       NFY_B_biallelic_peaks_CCAAT_shared,
                                                                       NFY_B_biallelic_peaks_jf1v2_shared,
                                                                       NFY_B_biallelic_peaks_CCAAT_jf1v2_shared,
                                                                       NFY_B_biallelic_peaks_CCAAT_jf1vindel_uniqB6,
                                                                       NFY_J_biallelic_peaks_CCAAT_shared,
                                                                       NFY_J_biallelic_peaks_jf1v2_shared,
                                                                       NFY_J_biallelic_peaks_CCAAT_jf1v2_shared,
                                                                       NFY_J_biallelic_peaks_CCAAT_jf1vindel_uniqJF1,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiuniq,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_InputBiallelicfilt_BJ1InputBiuniq,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq,
                                                                       Alsp_NFY_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq,
                                                                       NFY_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq,
                                                                       NFY_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq,
                                                                       NFY_Biallelicfilt_sharedB6JF1_peaks_NRF1_motif_SNP_Indelsuniq,
                                                                       NFY_Biallelicfilt_sharedB6JF1_peaks_pub_NFY_uniqBi,
                                                                       AlspNFYddsNormcounts_Biallelicfilt_BJ1InputBiSNPuniq)

colnames(matrix_AlspNFYddsNormcounts_Biallelicfilt_features) <- c("Chr","Start","End","PeakID","Score","IDX.1","IDY.1","B_CCAAT","IDX.2","IDY.2","B_SNPs","IDX.3","IDY.3","B_CCAAT_SNPs","IDX.4","IDY.4","B_CCAAT_Indel","IDX.5","IDY.5","J_CCAAT","IDX.6","IDY.6","J_SNPs","IDX.7","IDY.7","J_CCAAT_SNPs","IDX.8","IDY.8","J_CCAAT_Indel","IDX.9","IDY.9","sharedB6JF1_Indel","IDX.10","IDY.10","Input_Het_SNPs","IDX.11","IDY.11","Input_Bi_Regions","IDX.12","IDY.12","Input_Bi_Het_SNPs","IDX.13","IDY.13","ICRs","IDX.14","IDY.14","CGIs","IDX.15","IDY.15","SP1","IDX.16","IDY.16","YY1","IDX.17","IDY.17","NRF1","IDX.18","IDY.18","pub_uniqBi","IDX.19","PosHetSNPs")

dim(matrix_AlspNFYddsNormcounts_Biallelicfilt_features)
head(matrix_AlspNFYddsNormcounts_Biallelicfilt_features)
tail(matrix_AlspNFYddsNormcounts_Biallelicfilt_features)
matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab  <- data.frame(matrix_AlspNFYddsNormcounts_Biallelicfilt_features[,c(1:4,seq(5,59,3),61)])
matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab["Support_Index"] <- rowMeans(matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab[,6:(length(colnames(matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab))-1)])
matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab <- matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab[order(-matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab)

write.table(matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab, "matrix_AlspNFYddsNormcounts_Biallelicfilt_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)


#Upload reads to UCSC
#Upload reads for visulaization
sort -k1,1 -k2,2n /media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_B6_reads.sort.bed
sort -k1,1 -k2,2n  /media/ankitv/Archivio1/2021/chip_seq/BJ1_NFY/newway/allele_sp/Input/Input_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_JF1_reads.sort.bed


~/tools_av/bedToBigBed Input_B6_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_B6_reads.bb
~/tools_av/bedToBigBed Input_JF1_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_JF1_reads.bb

scp -r *.bb ankits@140.164.60.71:/media/ankits/Archivio_2/ankit/2021/chip-seq/Bigwigs/2021_new_chip-seq_peaks/
  
track type=bigBed name="BJ1_Input_reads" description="BJ1_Input_reads" bigDataUrl=http://140.164.60.71:8080/2021_new_chip-seq_peaks/Input_reads.bb visibility=dense color=0,0,0
track type=bigBed name="BJ1_Input_B6_reads" description="BJ1_Input_B6_reads" bigDataUrl=http://140.164.60.71:8080/2021_new_chip-seq_peaks/Input_B6_reads.bb visibility=dense color=0,0,0
track type=bigBed name="BJ1_Input_JF1_reads" description="BJ1_Input_JF1_reads" bigDataUrl=http://140.164.60.71:8080/2021_new_chip-seq_peaks/Input_JF1_reads.bb visibility=dense color=0,0,0


  
  
######################################################   END OF ANALYSIS  ##################################################à

  
  #Intersect BJ1 an BJ1
grep CCAAT -v  ./../BJ1_NFY/Matrix/matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab.txt | sort -k1,1 -k2,2n > BJ1_NFY_B6.txt
grep CCAAT -v  ./../BJ1_NFY/Matrix/matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab.txt | sort -k1,1 -k2,2n > BJ1_NFY_JF1.txt
grep CCAAT -v  ./../NFY_NRF1/way3_flash7/mm10/Matrix/matrix_AlspNFYddsNormcounts_B6_monofilt_uniqB6_features_tab.txt  | sort -k1,1 -k2,2n > BJ1_NFY_B6.txt
grep CCAAT -v  ./../NFY_NRF1/way3_flash7/mm10/Matrix/matrix_AlspNFYddsNormcounts_JF1_monofilt_uniqJF1_features_tab.txt  | sort -k1,1 -k2,2n > BJ1_NFY_JF1.txt

bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt |  wc -l
bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt |  wc -l
bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt |  wc -l
bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt |  wc -l
bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt | grep chrX -v | wc -l
bedtools intersect -wa -wb -a BJ1_NFY_B6.txt -b BJ1_NFY_B6.txt |  wc -l
bedtools intersect -wa -wb -a BJ1_NFY_JF1.txt -b BJ1_NFY_JF1.txt |  wc -l

bedtools intersect -wa -a BJ1_NFY_B6.txt -b BJ1_NFY_B6.txt >  NFY_BJ1_B6_BJ1_B6.txt
bedtools intersect -wa -a BJ1_NFY_JF1.txt -b BJ1_NFY_JF1.txt > NFY_BJ1_JF1_BJ1_JF1.txt
bedtools intersect -wa -a BJ1_NFY_B6.txt -b BJ1_NFY_JF1.txt >  NFY_BJ1_B6_BJ1_JF1.txt
bedtools intersect -wa -a BJ1_NFY_JF1.txt -b BJ1_NFY_B6.txt > NFY_BJ1_JF1_BJ1_B6.txt
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
bedtools coverage -a NFY_intersection_rep1_2_peaks.bed -b BJ1resultInputBi_chr.txt  > NFY_intersection_rep1_2_peaks_covSNP.bed
bedtools intersect -wa -wb -a NFY_intersection_rep1_2_peaks_covSNP.bed -b AlspNFYddsNormcounts_B6_monofilt_uniqB6.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""B6"}' > NFY_intersection_rep1_2_peaks_covSNP_B6.bed
bedtools intersect -wa -wb -a NFY_intersection_rep1_2_peaks_covSNP.bed -b AlspNFYddsNormcounts_JF1_monofilt_uniqJF1.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""JF1"}' > NFY_intersection_rep1_2_peaks_covSNP_JF1.bed
cat NFY_intersection_rep1_2_peaks_covSNP_B6.bed NFY_intersection_rep1_2_peaks_covSNP_JF1.bed > NFY_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed
bedtools intersect -wa -wb -a NFY_intersection_rep1_2_peaks_covSNP.bed -b NFY_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Biallelic"}' > NFY_intersection_rep1_2_peaks_covSNP_Biallelic.bed
cat NFY_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed NFY_intersection_rep1_2_peaks_covSNP_Biallelic.bed > NFY_intersection_rep1_2_peaks_covSNP_B6_or_JF1_or_Biallelic.bed

BJ1_NFY_peaks_Het_SNP <- read.table("NFY_intersection_rep1_2_peaks_covSNP_B6_or_JF1_or_Biallelic.bed", sep = "\t", header = F, stringsAsFactors = F)
BJ1_NFY_peaks_Het_SNP <- BJ1_NFY_peaks_Het_SNP[,c(1,2,3,6,7)]
head(BJ1_NFY_peaks_Het_SNP)
colnames(BJ1_NFY_peaks_Het_SNP) <- c("chr","start","end","SNP","allele")
head(BJ1_NFY_peaks_Het_SNP)
BJ1_NFY_peaks_Het_SNP["SNP_allele"] <- paste0(BJ1_NFY_peaks_Het_SNP$SNP, "_" ,BJ1_NFY_peaks_Het_SNP$allele)
countBJ1_NFY_peaks_Het_SNP <- count(BJ1_NFY_peaks_Het_SNP,"SNP_allele")
countBJ1_NFY_peaks_Het_SNPsplit <- cSplit(countBJ1_NFY_peaks_Het_SNP,"SNP_allele","_")
countBJ1_NFY_peaks_Het_SNPsplit <- countBJ1_NFY_peaks_Het_SNPsplit[order(countBJ1_NFY_peaks_Het_SNPsplit$SNP_allele_2),]

ggplot(data=countBJ1_NFY_peaks_Het_SNPsplit, aes(x=SNP_allele_1, y=freq)) +
  geom_bar(stat="identity", aes(fill=as.character(SNP_allele_2)),position=position_dodge())+
  theme_bw()+
  scale_fill_manual(values = c("darkred","darkgreen","navy"))
ggsave("countBJ1_NFY_peaks_Het_SNP.png", width=10*1.25, height=10*1.25, units="cm", dpi=96)

head(BJ1result[["BJ1"]])
SNPshetBJ1result <- BJ1result[["BJ1"]]
SNPshetBJ1result$REF.counts <- SNPshetBJ1result$REF.counts/2
SNPshetBJ1result$ALT.counts <- SNPshetBJ1result$ALT.counts/2

rownames(SNPshetBJ1result) <- paste0(SNPshetBJ1result$CHROM, "_", SNPshetBJ1result$POS)
SNPshetBJ1result <- SNPshetBJ1result[,c(6,7)]
SNPshetBJ1resultst <- stack(as.matrix(SNPshetBJ1result))
SNPshetBJ1resultst <- data.frame(SNPshetBJ1resultst)
ggplot(SNPshetBJ1resultst, aes(x=value, fill=col, color=col)) +
  geom_histogram(position="dodge",alpha=0.5)+theme_bw()
ggsave("SNPshetBJ1resultst.png", width=17*1.25, height=10*1.25, units="cm", dpi=96)
