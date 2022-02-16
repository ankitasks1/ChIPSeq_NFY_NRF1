print("Identification of Allele Specific Peaks")
print("Factor in Analysis :  NRF1")
print("Mouse Genome : version mm10")
print("Mouse cell line : JB1")
print("Raw data processing pipeline : ChIP-Seq")
library(DESeq2)
library(EDASeq)
library(metaseqR)
library(splitstackshape)
library(ggplot2)
library(plyr)
library(BaalChIP)
setwd("/media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/cutoff/enrichment_R1_NRF1")
#Which peak caller is good macs1 and macs2
bedtools intersect -a NRF1_rep1_peaks.bed -b ./../PeakCall_R1_MACS2/NRF1_rep1_peaks.narrowPeak | wc -l
10856
#Both are ok
#But I used macs1 output becuase it has slightly higher number of peaks
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NRF1/rep1/PeakCall_R1_MACS1/NRF1_rep1_peaks.bed ./
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NRF1/rep2/PeakCall_R1_MACS1/NRF1_rep2_peaks.bed ./
  
#Count Peaks
wc -l NRF1_rep1_peaks.bed 
wc -l NRF1_rep2_peaks.bed

#Index the bams which can be used in visualization
samtools index ./../../bulk/Input/Input_rep1_R1_uniq.rmdup.bam 
samtools index ./../../bulk/Input/Input_rep2_R1_uniq.rmdup.bam 
samtools index ./../../bulk/NRF1/rep1/NRF1_rep1_R1_uniq.rmdup.bam 
samtools index ./../../bulk/NRF1/rep2/NRF1_rep2_R1_uniq.rmdup.bam 
samtools index ./../../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam
samtools index ./../../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam
samtools index ./../../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam
samtools index ./../../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam
samtools index ./../../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam


#Bulk peak intersection (using only -a and -b option will aloow to intersect peaks in such a way that only shared region between two peaks will be printed)
bedtools intersect -a NRF1_rep1_peaks.bed -b NRF1_rep2_peaks.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"NR"\t"$5}' > NRF1_intersection_rep1_2_peaks.bed
#I selected No -u and let shared  intervals present and given NUMBER to each peak ID
#Coverage calculation, 3 possible ways
bedtools multicov -bams ./../../bulk/Input/Input_rep1_R1_uniq.rmdup.bam -bed NRF1_intersection_rep1_2_peaks.bed > Input_rep1_R1_NRF1_int_peaks.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../bulk/Input/Input_rep1_R1_uniq.rmdup.bam  > Input_rep1_NRF1_intersection_rep1_2_peaks.bed
#and featurecounts based
#For Bulk already contain chr in chromosme files

#But allele specific data does not contain chr so chr need to be added before peak calling and featurecounts
#First convert all picard dedup bam to bed (NOTE: When I convert picard output bed '.rmdup.' is used while manualy dedup file has '_rmdup.'. So this is the difference and can be seen in mnauscript)
#eg. Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed (picard)
#Input_rep1_R1_uniq.sortedByReadname.genome.sort_JF1_rmdup.bed

#Run BamtoBed in respective folders
#bedtools bamtobed -i Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
#bedtools bamtobed -i Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed

bedtools bamtobed -i NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
bedtools bamtobed -i NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed

#bedtools bamtobed -i Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
#bedtools bamtobed -i Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed

bedtools bamtobed -i NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed
bedtools bamtobed -i NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bed > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bed > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed


awk '{print $4"\t"$1"\t"$2"\t"$3"\t""."}' NRF1_intersection_rep1_2_peaks.bed > NRF1_intersection_rep1_2_peaks.saf
#/home/ankitv/tools_av/subread-2.0.0-source/bin/featureCounts -a NRF1_intersection_rep1_2_peaks.saf -F SAF -o JB1_NRF1_featurecounts_R1.txt ./../bulk/Input/Input_rep1_R1_uniq.rmdup.bam  ./../bulk/Input/Input_rep2_R1_uniq.rmdup.bam  ./../bulk/NRF1/rep1/NRF1_rep1_R1_uniq.rmdup.bam  ./../bulk/NRF1/rep2/NRF1_rep2_R1_uniq.rmdup.bam  ./../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam ./../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam ./../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam ./../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam ./../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam


NRF1_intersection_rep1_2_peaks <- read.table("NRF1_intersection_rep1_2_peaks.bed", header = F)
colnames(NRF1_intersection_rep1_2_peaks) <- c("chr","start","end","peakID","score")

bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../bulk/Input/Input_rep1_R1_uniq.rmdup.bam > Input_rep1_R1_uniq.rmdup.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../bulk/Input/Input_rep2_R1_uniq.rmdup.bam > Input_rep2_R1_uniq.rmdup.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../bulk/NRF1/rep1/NRF1_rep1_R1_uniq.rmdup.bam > NRF1_rep1_R1_uniq.rmdup.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../bulk/NRF1/rep2/NRF1_rep2_R1_uniq.rmdup.bam > NRF1_rep2_R1_uniq.rmdup.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/NRF1/rep1/NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b ./../../allele_sp/NRF1/rep2/NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed

   
paste Input_rep1_R1_uniq.rmdup.NRF1_cov.bed Input_rep2_R1_uniq.rmdup.NRF1_cov.bed NRF1_rep1_R1_uniq.rmdup.NRF1_cov.bed NRF1_rep2_R1_uniq.rmdup.NRF1_cov.bed Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.NRF1_cov.bed NRF1_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.NRF1_cov.bed > JB1_NRF1_covdata_R1.txt
#testfilefeat <- read.table("JB1_NRF1_featurecounts_R1_re.txt")
#colnames(testfilefeat) <- c("peakID","chr","start","end","strand","length","Input_1","Input_2","NRF1_1","NRF1_2","Input_B6_1","Input_JF1_1","Input_B6_2","Input_JF1_2","NRF1_B6_1","NRF1_JF1_1","NRF1_B6_2","NRF1_JF1_2")
#So I used bedtools coverage based read assignment

NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1 <- read.table("JB1_NRF1_covdata_R1.txt", header = F)
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1)
dim(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1)

#Test if paste parallely
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(seq(4,length(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1),9))])
tail(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(seq(4,length(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1),9))])
#So the data is pasted parallely
#Now prepare coverage file
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1 <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1[,c(1:5,seq(6,length(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1),9))]
colnames(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1) <- c("chr","start","end","peakID","score","Input_1","Input_2","NRF1_1","NRF1_2","Input_B6_1","Input_JF1_1","Input_B6_2","Input_JF1_2","NRF1_B6_1","NRF1_JF1_1","NRF1_B6_2","NRF1_JF1_2")
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
tail(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
#match with head(testfilefeat) tail((testfilefeat))

head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
tail(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
dim(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
rownames(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1) <- paste0(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1$chr,"%",
                                                                   NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1$start,"%",
                                                                   NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1$end,"%",
                                                                   NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1$peakID,"%",
                                                                   NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1$score)
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1 <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,c(6:17)]

head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
colSums(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr["peaks"] <- rownames(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- cSplit(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr, "peaks", "%") #library(splitstackshape)
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
dim(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr[,c(13:17, 1:12)]
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr)
write.table(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr,"NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Use Deseq2 for calculating scaling factor from bulk
#Input normalization not required as the peaks are called on background of input
head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
Bulk_NRF1_countdata <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,3:4]
head(Bulk_NRF1_countdata)

bulkNRF1coldata <- read.table("bulkNRF1coldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bulkNRF1coldata) <- colnames(Bulk_NRF1_countdata)
head(bulkNRF1coldata)
bulkNRF1coldata <- bulkNRF1coldata[,c("condition","replicate")]
bulkNRF1coldata$condition <- factor(bulkNRF1coldata$condition)
bulkNRF1coldata$replicate <- factor(bulkNRF1coldata$replicate)
all(rownames(bulkNRF1coldata) == colnames(Bulk_NRF1_countdata)) #should print TRUE
#Low counts filter (filter only rowsums for NRF1 allele specific)
#NRF1keep <- rowSums(Bulk_NRF1_countdata) >= 10
#Bulk_NRF1_countdatafilt <- Bulk_NRF1_countdata[NRF1keep,]
#dim(Bulk_NRF1_countdatafilt)

#Filter using proportion function of NOISeq:

Replicates <- factor(bulkNRF1coldata$replicate)
dim(Bulk_NRF1_countdata)
Bulk_NRF1_countdatafilt <- NOISeq::filtered.data(dataset=Bulk_NRF1_countdata, factor=Replicates, method=3, norm=FALSE)
head(Bulk_NRF1_countdatafilt)
plotRLE(as.matrix(Bulk_NRF1_countdatafilt), outline=F)
write.table(Bulk_NRF1_countdatafilt, "Bulk_NRF1_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)


bulkNRF1dds <- DESeqDataSetFromMatrix(countData =Bulk_NRF1_countdatafilt, colData = bulkNRF1coldata, design = ~ replicate)
bulkNRF1dds
#No filtering
bulkNRF1ddscounts <- counts(bulkNRF1dds, normalized=FALSE)
head(bulkNRF1ddscounts)
dim(bulkNRF1ddscounts)
colSums(bulkNRF1ddscounts)
#View filtered count matrix: View(counts(bulkNRF1dds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
bulkNRF1ddsNorm <- estimateSizeFactors(bulkNRF1dds)
sizeFactors(bulkNRF1ddsNorm)

head(NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1)
Alsp_NRF1_countdata <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,9:12]
head(Alsp_NRF1_countdata)
gcoldata <- read.table("NRF1gcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(gcoldata)
gcoldata[,1]
rownames(gcoldata)=gcoldata[,1]
rownames(gcoldata)
colnames(gcoldata)
gcoldata = gcoldata[,-1]
gcoldata <- data.frame(gcoldata)
head(gcoldata)
gcoldata <- gcoldata[,c("gcondition","greplicate","gallele")]
rownames(gcoldata) <- colnames(Alsp_NRF1_countdata)
gcoldata$gcondition <- factor(gcoldata$gcondition)
gcoldata$greplicate <- factor(gcoldata$greplicate)
gcoldata$gallele <- factor(gcoldata$gallele)

all(rownames(gcoldata) == colnames(Alsp_NRF1_countdata)) #should print TRUE

#Low counts filter (filter only rowsums for NRF1 allele specific)
#alNRF1keep <- rowSums(Alsp_NRF1_countdata[,6:9]) >= 10
#Alsp_NRF1_countdatafilt <- Alsp_NRF1_countdata[alNRF1keep,]
ASReplicates <- factor(gcoldata$gallele)
dim(Alsp_NRF1_countdata)
Alsp_NRF1_countdatafilt <- NOISeq::filtered.data(dataset=Alsp_NRF1_countdata, factor=ASReplicates, method=3, norm=FALSE)
head(Alsp_NRF1_countdatafilt)
plotRLE(as.matrix(Alsp_NRF1_countdatafilt), outline=F)
write.table(Alsp_NRF1_countdatafilt, "Alsp_NRF1_countdatafilt.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

dim(Alsp_NRF1_countdatafilt)
AlspNRF1dds <- DESeqDataSetFromMatrix(countData =Alsp_NRF1_countdatafilt, colData = gcoldata, design = ~ gallele)
AlspNRF1dds
#No filtering
AlspNRF1ddscounts <- counts(AlspNRF1dds, normalized=FALSE)
head(AlspNRF1ddscounts)
dim(AlspNRF1ddscounts)
colSums(AlspNRF1ddscounts)

#Normalize by DESeq2 size factor, set size factoe using sizeFactors(bulkNRF1ddsNorm)
AlspNRF1ddsNorm <- AlspNRF1dds
sizeFactors(AlspNRF1ddsNorm) <- c(1.3628987, 1.3628987, 0.7337302, 0.7337302)
sizeFactors(AlspNRF1ddsNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
AlspNRF1ddsNormcounts <- counts(AlspNRF1ddsNorm, normalized=TRUE)
head(AlspNRF1ddsNormcounts)
write.table(AlspNRF1ddsNormcounts, "AlspNRF1ddsNormcounts.txt", sep="\t", quote=F, append = F)
boxplot(AlspNRF1ddsNormcounts, ylim=c(0,100))
AlspNRF1ddsNormcounts <- data.frame(AlspNRF1ddsNormcounts)
AlspNRF1ddsNormcounts["matAllelicRatio1"] <- (AlspNRF1ddsNormcounts$NRF1_JF1_1 / (AlspNRF1ddsNormcounts$NRF1_B6_1 + AlspNRF1ddsNormcounts$NRF1_JF1_1))
AlspNRF1ddsNormcounts["matAllelicRatio2"] <- (AlspNRF1ddsNormcounts$NRF1_JF1_2 / (AlspNRF1ddsNormcounts$NRF1_B6_2 + AlspNRF1ddsNormcounts$NRF1_JF1_2))
head(AlspNRF1ddsNormcounts)
dim(AlspNRF1ddsNormcounts)

AlspNRF1ddsNormcounts_chr <- AlspNRF1ddsNormcounts
AlspNRF1ddsNormcounts_chr["peaks"] <- rownames(AlspNRF1ddsNormcounts_chr)
head(AlspNRF1ddsNormcounts_chr)
AlspNRF1ddsNormcounts_chr <- cSplit(AlspNRF1ddsNormcounts_chr, "peaks", "%") #library(splitstackshape)
head(AlspNRF1ddsNormcounts_chr)
AlspNRF1ddsNormcounts_chr <- AlspNRF1ddsNormcounts_chr[,c(7:11, 1:6)]
head(AlspNRF1ddsNormcounts_chr)
colnames(AlspNRF1ddsNormcounts_chr) <- c("chr", "start", "end", "peakID", "score","NRF1_B6_1","NRF1_JF1_1","NRF1_B6_2","NRF1_JF1_2","matAllelicRatio1","matAllelicRatio2")
write.table(AlspNRF1ddsNormcounts_chr,"AlspNRF1ddsNormcounts_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

##Predict Monoallelic and Biallelic Peaks  --------------###
#------------ I used Proportion test -------------#
#Combined p.vales -2Sumation(log(Pi))
#ChiSquare df = 2*length(p-values), so if 2 p.values 2 *2
#lower.tail =FALSE and lower.tail =TRUE, 1-p
#Prop test and tag chromsome positions
#AlspNRF1ddsNormcounts_chr <- read.table("AlspNRF1ddsNormcounts_chr.txt", header = F)

head(AlspNRF1ddsNormcounts_chr)
dim(AlspNRF1ddsNormcounts_chr)
AlspNRF1ddsNormcounts_chr <- data.frame(AlspNRF1ddsNormcounts_chr)
rownames(AlspNRF1ddsNormcounts_chr) <- paste(AlspNRF1ddsNormcounts_chr$chr,AlspNRF1ddsNormcounts_chr$start,AlspNRF1ddsNormcounts_chr$end,AlspNRF1ddsNormcounts_chr$peakID,AlspNRF1ddsNormcounts_chr$score, sep = "%")
AlspNRF1ddsNormcounts_chr_prop <- AlspNRF1ddsNormcounts_chr[,6:11]
head(AlspNRF1ddsNormcounts_chr_prop)
#Take sum of alleles
AlspNRF1ddsNormcounts_chr_prop["NRF1_alavg_1"] <- AlspNRF1ddsNormcounts_chr_prop$NRF1_B6_1 + AlspNRF1ddsNormcounts_chr_prop$NRF1_JF1_1
AlspNRF1ddsNormcounts_chr_prop["NRF1_alavg_2"] <- AlspNRF1ddsNormcounts_chr_prop$NRF1_B6_2 + AlspNRF1ddsNormcounts_chr_prop$NRF1_JF1_2

#Take allelic ratio
AlspNRF1ddsNormcounts_chr_prop["mat_allele_ratio"] <- (AlspNRF1ddsNormcounts_chr_prop$matAllelicRatio1 + AlspNRF1ddsNormcounts_chr_prop$matAllelicRatio2)/2


dim(AlspNRF1ddsNormcounts_chr_prop) #2762    9
AlspNRF1ddsNormcounts_chr_propid <- AlspNRF1ddsNormcounts_chr_prop
AlspNRF1ddsNormcounts_chr_propid["id"] <- rownames(AlspNRF1ddsNormcounts_chr_prop)
head(AlspNRF1ddsNormcounts_chr_propid)
writexl::write_xlsx(AlspNRF1ddsNormcounts_chr_propid, "AlspNRF1ddsNormcounts_chr_prop.xlsx")

dim(AlspNRF1ddsNormcounts_chr_prop)#2762    9
head(AlspNRF1ddsNormcounts_chr_prop,1)
summary(AlspNRF1ddsNormcounts_chr_prop)
AlspNRF1ddsNormcounts_chr_prop <- data.frame(AlspNRF1ddsNormcounts_chr_prop)
#NA detected as some in some cases both allele were 0, So That need to be removed
AlspNRF1ddsNormcounts_chr_prop.filt <- na.omit(AlspNRF1ddsNormcounts_chr_prop) 
head(AlspNRF1ddsNormcounts_chr_prop.filt)
dim(AlspNRF1ddsNormcounts_chr_prop.filt)
summary(AlspNRF1ddsNormcounts_chr_prop.filt)
NRF1_alavg_1.prop1 <- Map(prop.test,x =AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_B6_1, n= AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_1, p=0.5)

AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg_1.pvalue"] <- data.frame(capture.output(for (i in seq_along(NRF1_alavg_1.prop1)){
  cat(NRF1_alavg_1.prop1[[i]]$p.value, "\n")
}))

NRF1_alavg_2.prop1 <- Map(prop.test,x =AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_B6_2, n= AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_2, p=0.5)

AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg_2.pvalue"] <- data.frame(capture.output(for (i in seq_along(NRF1_alavg_2.prop1)){
  cat(NRF1_alavg_2.prop1[[i]]$p.value, "\n")
}))


#Combine p-values calculated from prop.test for rep1 and rep2
AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg.comb.testpvalue"] <- -2 *(log(as.numeric(as.character(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_1.pvalue))) + log(as.numeric(as.character(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_2.pvalue))))
#To deal with -2 * (log(0)+log(0)), I used fisher.method function which has option to zero.sub to value, I set it 2.2e-16, this will replace Inf with 2.2e-16, 
#eg. -2*(log(2.2e-16)+log(2.2e-16)) = 144.2116
class(AlspNRF1ddsNormcounts_chr_prop.filt)
AlspNRF1ddsNormcounts_chr_prop.filt <- data.frame(AlspNRF1ddsNormcounts_chr_prop.filt)
AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg.comb.fisherpvalue"] <- fisher.method(data.frame(as.numeric(as.character(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_1.pvalue)), as.numeric(as.character(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg_2.pvalue))), method = c("fisher"), p.corr ="none", zero.sub = 2.2e-16, na.rm = FALSE, mc.cores=NULL)
#Warning will pop up Warning message:
#In `[<-.data.frame`(`*tmp*`, "NRF1_alavg.comb.fisherpvalue", value = list( :provided 4 variables to replace 1 variables
#Export the sheet and manually checked the column NRF1_alavg.comb.testpvalue and NRF1_alavg.comb.fisherpvalue should match except Inf, so manual and fisher.method is ok
head(AlspNRF1ddsNormcounts_chr_prop.filt)
AlspNRF1ddsNormcounts_chr_prop.filt_check <- AlspNRF1ddsNormcounts_chr_prop.filt
AlspNRF1ddsNormcounts_chr_prop.filt_check["ID"] <- rownames(AlspNRF1ddsNormcounts_chr_prop.filt_check)
writexl::write_xlsx(AlspNRF1ddsNormcounts_chr_prop.filt_check, "AlspNRF1ddsNormcounts_chr_prop.filt_check.xlsx")

head(AlspNRF1ddsNormcounts_chr_prop.filt)
tail(AlspNRF1ddsNormcounts_chr_prop.filt)
#Chisquare test for p-values#pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)#p vector of probabilities, n	 number of observations, degrees of freedom (non-negative, but can be non-integer)
AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg.comb.pvalue.pchisq"] <- pchisq(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg.comb.fisherpvalue,4, lower.tail=FALSE)
AlspNRF1ddsNormcounts_chr_prop.filt["NRF1_alavg.comb.fdr"] <- p.adjust(AlspNRF1ddsNormcounts_chr_prop.filt$NRF1_alavg.comb.pvalue.pchisq,method="BH")
head(AlspNRF1ddsNormcounts_chr_prop.filt)

#Introduce peak info
AlspNRF1ddsNormcounts_chr_prop.filt["id"] <- data.frame(rownames(AlspNRF1ddsNormcounts_chr_prop.filt))
head(AlspNRF1ddsNormcounts_chr_prop.filt)
AlspNRF1ddsNormcounts_chr_prop.filt.sep <- cSplit(AlspNRF1ddsNormcounts_chr_prop.filt, "id", "%")
dim(AlspNRF1ddsNormcounts_chr_prop.filt.sep)
head(AlspNRF1ddsNormcounts_chr_prop.filt.sep)
AlspNRF1ddsNormcounts_chr_prop.filt.sep <- AlspNRF1ddsNormcounts_chr_prop.filt.sep[,c(16:20,1:15)]
colnames(AlspNRF1ddsNormcounts_chr_prop.filt.sep) <- c("chr","start","end","peakID","score",colnames(AlspNRF1ddsNormcounts_chr_prop.filt[,1:15]))
#JF1_monoallelic
AlspNRF1ddsNormcounts_JF1_mono <- AlspNRF1ddsNormcounts_chr_prop.filt.sep[which(AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.75 &
                                                                                AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.75),]
head(AlspNRF1ddsNormcounts_JF1_mono)
dim(AlspNRF1ddsNormcounts_JF1_mono)


#Monoalleically expressed to JF1 with FDR < 0.05
AlspNRF1ddsNormcounts_JF1_monofilt <- AlspNRF1ddsNormcounts_JF1_mono[which(AlspNRF1ddsNormcounts_JF1_mono$NRF1_alavg.comb.fdr<0.05),]
head(AlspNRF1ddsNormcounts_JF1_monofilt)
dim(AlspNRF1ddsNormcounts_JF1_monofilt)
write.table(AlspNRF1ddsNormcounts_JF1_monofilt,"AlspNRF1ddsNormcounts_JF1_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNRF1ddsNormcounts_JF1_monofilt_qual <- AlspNRF1ddsNormcounts_JF1_monofilt
head(AlspNRF1ddsNormcounts_JF1_monofilt_qual)
AlspNRF1ddsNormcounts_JF1_monofilt_qual199 <- AlspNRF1ddsNormcounts_JF1_monofilt_qual[which(AlspNRF1ddsNormcounts_JF1_monofilt_qual$score >= 199),]
dim(AlspNRF1ddsNormcounts_JF1_monofilt_qual199)
head(AlspNRF1ddsNormcounts_JF1_monofilt_qual199)
write.table(AlspNRF1ddsNormcounts_JF1_monofilt_qual199,"AlspNRF1ddsNormcounts_JF1_monofilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#B6_monoallelic
AlspNRF1ddsNormcounts_B6_mono <- AlspNRF1ddsNormcounts_chr_prop.filt.sep[which(AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.25 &
                                                                               AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.25),]
head(AlspNRF1ddsNormcounts_B6_mono)
dim(AlspNRF1ddsNormcounts_B6_mono)

#Monoalleically expressed to B6 with FDR < 0.05
AlspNRF1ddsNormcounts_B6_monofilt <- AlspNRF1ddsNormcounts_B6_mono[which(AlspNRF1ddsNormcounts_B6_mono$NRF1_alavg.comb.fdr<0.05),]
head(AlspNRF1ddsNormcounts_B6_monofilt)
dim(AlspNRF1ddsNormcounts_B6_monofilt)
write.table(AlspNRF1ddsNormcounts_B6_monofilt,"AlspNRF1ddsNormcounts_B6_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNRF1ddsNormcounts_B6_monofilt_qual <- AlspNRF1ddsNormcounts_B6_monofilt
head(AlspNRF1ddsNormcounts_B6_monofilt_qual)
AlspNRF1ddsNormcounts_B6_monofilt_qual199 <- AlspNRF1ddsNormcounts_B6_monofilt_qual[which(AlspNRF1ddsNormcounts_B6_monofilt_qual$score >= 199),]
dim(AlspNRF1ddsNormcounts_B6_monofilt_qual199)
head(AlspNRF1ddsNormcounts_B6_monofilt_qual199)
write.table(AlspNRF1ddsNormcounts_B6_monofilt_qual199,"AlspNRF1ddsNormcounts_B6_monofilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#Biallelic
#Biallelicalleleic
AlspNRF1ddsNormcounts_Biallelic <- AlspNRF1ddsNormcounts_chr_prop.filt.sep[which(AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.60 &
                                                                                 AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.60 &
                                                                                 AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.40 &
                                                                                 AlspNRF1ddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.40),]
head(AlspNRF1ddsNormcounts_Biallelic)
dim(AlspNRF1ddsNormcounts_Biallelic)

#Bialleically expressed  with FDR > 0.05
AlspNRF1ddsNormcounts_Biallelicfilt <- AlspNRF1ddsNormcounts_Biallelic[which(AlspNRF1ddsNormcounts_Biallelic$NRF1_alavg.comb.fdr>0.05),]
head(AlspNRF1ddsNormcounts_Biallelicfilt)
dim(AlspNRF1ddsNormcounts_Biallelicfilt)
write.table(AlspNRF1ddsNormcounts_Biallelicfilt,"AlspNRF1ddsNormcounts_Biallelicfilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Quality filter
AlspNRF1ddsNormcounts_Biallelicfilt_qual <- AlspNRF1ddsNormcounts_Biallelicfilt
head(AlspNRF1ddsNormcounts_Biallelicfilt_qual)
AlspNRF1ddsNormcounts_Biallelicfilt_qual199 <- AlspNRF1ddsNormcounts_Biallelicfilt_qual[which(AlspNRF1ddsNormcounts_Biallelicfilt_qual$score >= 199),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_qual199)
head(AlspNRF1ddsNormcounts_Biallelicfilt_qual199)
write.table(AlspNRF1ddsNormcounts_Biallelicfilt_qual199,"AlspNRF1ddsNormcounts_Biallelicfilt_qual199.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

##No need to overlap on allele specific called peaks
##bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt.txt -b ./../merge_peakanalysis/NRF1_alsp_uniqB6_bulk.bed | sort -k4,4 -u > AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt
##bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt.txt -b ./../merge_peakanalysis/NRF1_alsp_uniqJF1_bulk.bed | sort -k4,4 -u > AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt
##bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt.txt -b ./../merge_peakanalysis/NRF1_alsp_sharedB6JF1_bulk.bed | sort -k4,4 -u > AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt

cp AlspNRF1ddsNormcounts_B6_monofilt.txt AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt
cp AlspNRF1ddsNormcounts_JF1_monofilt.txt AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt
cp AlspNRF1ddsNormcounts_Biallelicfilt.txt AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt

##Overlap with ICR mm10 (for checking monoallelecity) 
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.sorted.bed > AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt


#Perform similar analysis on Input
#Purpose: Heterozygosity map using Input of Bulk and Allele specific
### ----------------Heterozygous Peaks in Input---------------------- ###

Bulk_Input_countdata <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,1:2]
head(Bulk_Input_countdata)
Alsp_Input_countdata <- NRF1_and_Input_intersection_R1_bulk_AS_B6_JF1_1[,5:8]
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
sizeFactors(AlspInputddsNorm) <- c(0.7405316, 0.7405316, 1.3503812, 1.3503812)
sizeFactors(AlspInputddsNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
AlspInputddsNormcounts <- counts(AlspInputddsNorm, normalized=TRUE)
head(AlspInputddsNormcounts)
write.table(AlspInputddsNormcounts, "AlspInputddsNormcounts.txt", sep="\t", quote=F, append = F)
boxplot(AlspInputddsNormcounts, ylim=c(0,100))
AlspInputddsNormcounts <- data.frame(AlspInputddsNormcounts)
AlspInputddsNormcounts["matAllelicRatio1"] <- (AlspInputddsNormcounts$Input_JF1_1 / (AlspInputddsNormcounts$Input_B6_1 + AlspInputddsNormcounts$Input_JF1_1))
AlspInputddsNormcounts["matAllelicRatio2"] <- (AlspInputddsNormcounts$Input_JF1_2 / (AlspInputddsNormcounts$Input_B6_2 + AlspInputddsNormcounts$Input_JF1_2))
head(AlspInputddsNormcounts)
dim(AlspInputddsNormcounts)

AlspInputddsNormcounts_chr <- AlspInputddsNormcounts
AlspInputddsNormcounts_chr["peaks"] <- rownames(AlspInputddsNormcounts_chr)
head(AlspInputddsNormcounts_chr)
AlspInputddsNormcounts_chr <- cSplit(AlspInputddsNormcounts_chr, "peaks", "%") #library(splitstackshape)
head(AlspInputddsNormcounts_chr)
AlspInputddsNormcounts_chr <- AlspInputddsNormcounts_chr[,c(7:11, 1:6)]
head(AlspInputddsNormcounts_chr)
colnames(AlspInputddsNormcounts_chr) <- c("chr", "start", "end", "peakID", "score","Input_B6_1","Input_JF1_1","Input_B6_2","Input_JF1_2","matAllelicRatio1","matAllelicRatio2")
#write.table(AlspInputddsNormcounts_chr,"AlspInputddsNormcounts_chr.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

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
AlspInputddsNormcounts_chr_prop <- AlspInputddsNormcounts_chr[,6:11]
head(AlspInputddsNormcounts_chr_prop)
#Take sum of alleles
AlspInputddsNormcounts_chr_prop["Input_alavg_1"] <- AlspInputddsNormcounts_chr_prop$Input_B6_1 + AlspInputddsNormcounts_chr_prop$Input_JF1_1
AlspInputddsNormcounts_chr_prop["Input_alavg_2"] <- AlspInputddsNormcounts_chr_prop$Input_B6_2 + AlspInputddsNormcounts_chr_prop$Input_JF1_2

#Take allelic ratio
AlspInputddsNormcounts_chr_prop["mat_allele_ratio"] <- (AlspInputddsNormcounts_chr_prop$matAllelicRatio1 + AlspInputddsNormcounts_chr_prop$matAllelicRatio2)/2


dim(AlspInputddsNormcounts_chr_prop) #2773    9
AlspInputddsNormcounts_chr_propid <- AlspInputddsNormcounts_chr_prop
AlspInputddsNormcounts_chr_propid["id"] <- rownames(AlspInputddsNormcounts_chr_prop)
head(AlspInputddsNormcounts_chr_propid)
#writexl::write_xlsx(AlspInputddsNormcounts_chr_propid, "AlspInputddsNormcounts_chr_prop.xlsx")

dim(AlspInputddsNormcounts_chr_prop)#2773    9
head(AlspInputddsNormcounts_chr_prop,1)
summary(AlspInputddsNormcounts_chr_prop)
AlspInputddsNormcounts_chr_prop <- data.frame(AlspInputddsNormcounts_chr_prop)
#NA detected as some in some cases both allele were 0, So That need to be removed
AlspInputddsNormcounts_chr_prop.filt <- na.omit(AlspInputddsNormcounts_chr_prop) 
head(AlspInputddsNormcounts_chr_prop.filt)
dim(AlspInputddsNormcounts_chr_prop.filt)
summary(AlspInputddsNormcounts_chr_prop.filt)
Input_alavg_1.prop1 <- Map(prop.test,x =AlspInputddsNormcounts_chr_prop.filt$Input_B6_1, n= AlspInputddsNormcounts_chr_prop.filt$Input_alavg_1, p=0.5)

AlspInputddsNormcounts_chr_prop.filt["Input_alavg_1.pvalue"] <- data.frame(capture.output(for (i in seq_along(Input_alavg_1.prop1)){
  cat(Input_alavg_1.prop1[[i]]$p.value, "\n")
}))

Input_alavg_2.prop1 <- Map(prop.test,x =AlspInputddsNormcounts_chr_prop.filt$Input_B6_2, n= AlspInputddsNormcounts_chr_prop.filt$Input_alavg_2, p=0.5)

AlspInputddsNormcounts_chr_prop.filt["Input_alavg_2.pvalue"] <- data.frame(capture.output(for (i in seq_along(Input_alavg_2.prop1)){
  cat(Input_alavg_2.prop1[[i]]$p.value, "\n")
}))


#Combine p-values calculated from prop.test for rep1 and rep2
AlspInputddsNormcounts_chr_prop.filt["Input_alavg.comb.testpvalue"] <- -2 *(log(as.numeric(as.character(AlspInputddsNormcounts_chr_prop.filt$Input_alavg_1.pvalue))) + log(as.numeric(as.character(AlspInputddsNormcounts_chr_prop.filt$Input_alavg_2.pvalue))))
#To deal with -2 * (log(0)+log(0)), I used fisher.method function which has option to zero.sub to value, I set it 2.2e-16, this will replace Inf with 2.2e-16, 
#eg. -2*(log(2.2e-16)+log(2.2e-16)) = 144.2116
class(AlspInputddsNormcounts_chr_prop.filt)
AlspInputddsNormcounts_chr_prop.filt <- data.frame(AlspInputddsNormcounts_chr_prop.filt)
AlspInputddsNormcounts_chr_prop.filt["Input_alavg.comb.fisherpvalue"] <- fisher.method(data.frame(as.numeric(as.character(AlspInputddsNormcounts_chr_prop.filt$Input_alavg_1.pvalue)), as.numeric(as.character(AlspInputddsNormcounts_chr_prop.filt$Input_alavg_2.pvalue))), method = c("fisher"), p.corr ="none", zero.sub = 2.2e-16, na.rm = FALSE, mc.cores=NULL)
#Warning will pop up Warning message:
#In `[<-.data.frame`(`*tmp*`, "Input_alavg.comb.fisherpvalue", value = list( :provided 4 variables to replace 1 variables
#Export the sheet and manually checked the column Input_alavg.comb.testpvalue and Input_alavg.comb.fisherpvalue should match except Inf, so manual and fisher.method is ok
head(AlspInputddsNormcounts_chr_prop.filt)
AlspInputddsNormcounts_chr_prop.filt_check <- AlspInputddsNormcounts_chr_prop.filt
AlspInputddsNormcounts_chr_prop.filt_check["ID"] <- rownames(AlspInputddsNormcounts_chr_prop.filt_check)
#writexl::write_xlsx(AlspInputddsNormcounts_chr_prop.filt_check, "AlspInputddsNormcounts_chr_prop.filt_check.xlsx")

head(AlspInputddsNormcounts_chr_prop.filt)
tail(AlspInputddsNormcounts_chr_prop.filt)
#Chisquare test for p-values#pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)#p vector of probabilities, n	 number of observations, degrees of freedom (non-negative, but can be non-integer)
AlspInputddsNormcounts_chr_prop.filt["Input_alavg.comb.pvalue.pchisq"] <- pchisq(AlspInputddsNormcounts_chr_prop.filt$Input_alavg.comb.fisherpvalue,4, lower.tail=FALSE)
AlspInputddsNormcounts_chr_prop.filt["Input_alavg.comb.fdr"] <- p.adjust(AlspInputddsNormcounts_chr_prop.filt$Input_alavg.comb.pvalue.pchisq,method="BH")
head(AlspInputddsNormcounts_chr_prop.filt)

#Introduce gene name
AlspInputddsNormcounts_chr_prop.filt["id"] <- data.frame(rownames(AlspInputddsNormcounts_chr_prop.filt))
head(AlspInputddsNormcounts_chr_prop.filt)
AlspInputddsNormcounts_chr_prop.filt.sep <- cSplit(AlspInputddsNormcounts_chr_prop.filt, "id", "%")
dim(AlspInputddsNormcounts_chr_prop.filt.sep)
head(AlspInputddsNormcounts_chr_prop.filt.sep)
AlspInputddsNormcounts_chr_prop.filt.sep <- AlspInputddsNormcounts_chr_prop.filt.sep[,c(16:20,1:15)]
colnames(AlspInputddsNormcounts_chr_prop.filt.sep) <- c("chr","start","end","peakID","score",colnames(AlspInputddsNormcounts_chr_prop.filt[,1:15]))
#JF1_monoallelic
AlspInputddsNormcounts_JF1_mono <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.75 &
                                                                                    AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.75),]
head(AlspInputddsNormcounts_JF1_mono)
dim(AlspInputddsNormcounts_JF1_mono)


#Monoalleically expressed to JF1 with FDR < 0.05
AlspInputddsNormcounts_JF1_monofilt <- AlspInputddsNormcounts_JF1_mono[which(AlspInputddsNormcounts_JF1_mono$Input_alavg.comb.fdr<0.05),]
head(AlspInputddsNormcounts_JF1_monofilt)
dim(AlspInputddsNormcounts_JF1_monofilt)
write.table(AlspInputddsNormcounts_JF1_monofilt,"AlspInputddsNormcounts_JF1_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)

#B6_monoallelic
AlspInputddsNormcounts_B6_mono <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.25 &
                                                                                   AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.25),]
head(AlspInputddsNormcounts_B6_mono)
dim(AlspInputddsNormcounts_B6_mono)

#Monoalleically expressed to B6 with FDR < 0.05
AlspInputddsNormcounts_B6_monofilt <- AlspInputddsNormcounts_B6_mono[which(AlspInputddsNormcounts_B6_mono$Input_alavg.comb.fdr<0.05),]
head(AlspInputddsNormcounts_B6_monofilt)
dim(AlspInputddsNormcounts_B6_monofilt)
write.table(AlspInputddsNormcounts_B6_monofilt,"AlspInputddsNormcounts_B6_monofilt.txt", sep="\t",row.names=FALSE,col.names = FALSE, quote=FALSE, append = FALSE)


#Biallelic
#Biallelicalleleic
AlspInputddsNormcounts_Biallelic <- AlspInputddsNormcounts_chr_prop.filt.sep[which(AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 <= 0.60 &
                                                                                     AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 <= 0.60 &
                                                                                     AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio1 >= 0.40 &
                                                                                     AlspInputddsNormcounts_chr_prop.filt.sep$matAllelicRatio2 >= 0.40),]
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
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/mouse/motif[GCGC].mm10.bed GCGC_motif.mm10.sort.bed
#GCGC and JF1-B6 SNP
bedtools intersect -wa -wb -a GCGC_motif.mm10.sort.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > GCGC_jf1v2SNPmm10.bed

#Motif and SNP overlap peaks
#GCGC
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b GCGC_motif.mm10.sort.bed > NRF1_B6_specific_peaks_GCGC.bed

#SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_B6_specific_peaks_jf1v2.bed

#Note GCGC_jf1v2SNPmm10.bed is file intersected with SNP (jf1v2mm10)
#GCGC+SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b GCGC_jf1v2SNPmm10.bed > NRF1_B6_specific_peaks_GCGC_jf1v2.bed

#GCGC, SNP
bedtools intersect -wa -wb -a NRF1_B6_specific_peaks_GCGC.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_B6_specific_peaks_GCGC_then_jf1v2.bed

#NO NRF1
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b GCGC_motif.mm10.sort.bed -v > Nomotif_NRF1_B6_specific_peaks.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b GCGC_jf1v2SNPmm10.bed -v > NomotifSNP_NRF1_B6_specific_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NRF1_B6_specific_peaks.bed

#Monoallelic JF1
#Create a motif overlapped file JF1
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/mouse/motif[GCGC].mm10.jf1.bed GCGC_motif.mm10.jf1.sort.bed
#GCGC and JF1-JF1 SNP
bedtools intersect -wa -wb -a GCGC_motif.mm10.jf1.sort.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > GCGC_jf1v2SNPmm10.jf1.bed

#Motif and SNP overlap peaks
#GCGC
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b GCGC_motif.mm10.jf1.sort.bed > NRF1_JF1_specific_peaks_GCGC.bed

#SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_JF1_specific_peaks_jf1v2.bed

#Note GCGC_jf1v2SNPmm10.jf1.bed is file intersected with SNP (jf1v2mm10.jf1)
#GCGC+SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b GCGC_jf1v2SNPmm10.jf1.bed > NRF1_JF1_specific_peaks_GCGC_jf1v2.bed

#GCGC, SNP
bedtools intersect -wa -wb -a NRF1_JF1_specific_peaks_GCGC.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_JF1_specific_peaks_GCGC_then_jf1v2.bed

#NO NRF1
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b GCGC_motif.mm10.jf1.sort.bed -v > Nomotif_NRF1_JF1_specific_peaks.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b GCGC_jf1v2SNPmm10.jf1.bed -v > NomotifSNP_NRF1_JF1_specific_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NRF1_JF1_specific_peaks.bed

#Indel
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndel.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel.txt

#GCGC and JF1-B6 Indels
bedtools intersect -wa -wb -a GCGC_motif.mm10.sort.bed -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed > GCGC_jf1vindelSNPmm10.bed
bedtools intersect -wa -wb -a GCGC_motif.mm10.jf1.sort.bed -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed > GCGC_jf1vindelSNPmm10.jf1.bed

#GCGC+Indels
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b GCGC_jf1vindelSNPmm10.bed > NRF1_B6_specific_peaks_GCGC_jf1vindel.bed

#GCGC+Indels
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b GCGC_jf1vindelSNPmm10.jf1.bed > NRF1_JF1_specific_peaks_GCGC_jf1vindel.bed

#4321634 /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed
#32233 GCGC_jf1vindelSNPmm10.bed
#32045 GCGC_jf1vindelSNPmm10.jf1.bed



#------------------------- Identify Heterozygous SNPs  ------------------------#
#Allelic Bias count was obtained using BaalChIP
#Get Peak coordinates summit +/-100bp
filepath <- "/media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/cutoff/enrichment_R1_NRF1"
awk '{print $1,$2-100,$3+100,$4,$5}' OFS='\t' NRF1_rep1_summits.bed > NRF1_rep1_summits_span100.bed
awk '{print $1,$2-100,$3+100,$4,$5}' OFS='\t' NRF1_rep2_summits.bed > NRF1_rep2_summits_span100.bed

#Share span extended summits in rep1 and rep2
bedtools intersect -a NRF1_rep1_summits_span100.bed -b NRF1_rep2_summits_span100.bed  > NRF1_intersection_rep1_2_summits_span100.bed 

#Get SNPs overlapped summit +/- 100bp
bedtools intersect -wa -wb -a NRF1_intersection_rep1_2_summits_span100.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed | awk '{print $6,$7,$8,$9}' OFS='_' | sort -k1,1 -u | awk -F'_' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}'  > jf1v2_Snp_GRm10_NRF1summitsp100.bed
awk '{print NR"\t"$1"\t"$2"\t"$4"\t"$5}' jf1v2_Snp_GRm10_NRF1summitsp100.bed > jf1v2_hetSNP.txt
vim jf1v2_hetSNP.txt
ID	CHROM	POS	REF	ALT
jf1v2hets <- c("JB1"="jf1v2_hetSNP.txt")
head(read.delim(file.path(filepath,"jf1v2_hetSNP.txt")))
#Copy summit span to other name
cp NRF1_intersection_rep1_2_summits_span100.bed JB1_NRF1.bed
mysamplesheet <- file.path(filepath,"baalchipsamplesheet.tsv")

#mysamplesheet <- read.delim(file.path(filepath,"mysamplesheet.tsv"))
mysamplesheet

JB1res <- BaalChIP(samplesheet=mysamplesheet, hets=jf1v2hets)
JB1res
BaalChIP.get(JB1res, what="samples")
#run alleleCounts
JB1res <- alleleCounts(JB1res, min_base_quality=10, min_mapq=15, verbose=FALSE)

#run QC filter
#res <- QCfilter(res, 
#                RegionsToFilter=list("blacklist"=blacklist_mm10), 
#                verbose=FALSE)
JB1res <- mergePerGroup(JB1res)
JB1res <- filter1allele(JB1res)
JB1res <- getASB(JB1res, Iter=5000, conf_level=0.95, cores = 8, 
                 RMcorrection = FALSE, 
                 RAFcorrection=FALSE)
Sys.time()
JB1result <- BaalChIP.report(JB1res)
head(JB1result[["JB1"]])
#show ASB SNPs
JB1resultInput <- data.frame(JB1result[["JB1"]])
head(JB1resultInput)
dim(JB1resultInput)
JB1resultInput <- JB1resultInput[which(JB1resultInput$REF.counts >= 2),]
JB1resultInput <- JB1resultInput[which(JB1resultInput$ALT.counts >= 2),]

write.table(JB1resultInput, "JB1resultInputBi.txt", sep = "\t", append = F, quote = F, row.names = F)

JB1resultInputASBSNP <- JB1resultInput[JB1resultInput$isASB==TRUE,]
head(JB1resultInputASBSNP)
dim(JB1resultInputASBSNP)
write.table(JB1resultInputASBSNP, "JB1resultInputASBSNP.txt", sep = "\t", append = F, quote = F, row.names = F)

JB1resultInputBi <- read.table("JB1resultInputBi.txt", header = T, stringsAsFactors = F)
head(JB1resultInputBi)
dim(JB1resultInputBi)
JB1resultInputBi_chr <- JB1resultInputBi[,c(2,3,3:15)]
head(JB1resultInputBi_chr)
write.table(JB1resultInputBi_chr,  "JB1resultInputBi_chr.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)


#Heterozygous SNPs file: JB1resultInputBi_chr.txt
#Heterozygouse/Biallleic region in Input  file: AlspInputddsNormcounts_Biallelicfilt.txt
#Add biases peak in input information to matrix 

bedtools intersect -wa -wb -a AlspInputddsNormcounts_Biallelicfilt.txt -b JB1resultInputBi_chr.txt > AlspInputddsNormcounts_Biallelicfilt_JB1resultInputBi.txt

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b JB1resultInputBi_chr.txt > AlspNRF1ddsNormcounts_B6_JB1InputBi.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b JB1resultInputBi_chr.txt > AlspNRF1ddsNormcounts_JF1_JB1InputBi.txt

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNRF1ddsNormcounts_B6_InputBiallelicfilt.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt.txt


bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b AlspInputddsNormcounts_Biallelicfilt_JB1resultInputBi.txt > AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b AlspInputddsNormcounts_Biallelicfilt_JB1resultInputBi.txt > AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi.txt


#------------ The further characteristics of peaks are related to regions +/- 100 from the middle of peaks
#Peaks at +/- 100bp from center of peaks
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6 <- read.table("AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100 <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6
head(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100)
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100["diff"] <- round((AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$V3-AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$V2)/2)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100["chr"] <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$V1
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100["start"] <- (AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$V2 + AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$diff)-100
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100["end"] <- (AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$V3 - AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100$diff)+100
head(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100 <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100[,c(22:24,4,5)]
write.table(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100, "AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)

AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1 <- read.table("AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100 <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1
head(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100["diff"] <- round((AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V3-AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V2)/2)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100["chr"] <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V1
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100["start"] <- (AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V2 + AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$diff)-100
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100["end"] <- (AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$V3 - AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100$diff)+100
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100 <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100[,c(22:24,4,5)]
write.table(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100, "AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)

AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1 <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100 <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1
head(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["diff"] <- round((AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V3-AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V2)/2)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["chr"] <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V1
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["start"] <- (AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V2 + AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$diff)-100
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100["end"] <- (AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$V3 - AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100$diff)+100
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100 <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100[,c(22:24,4,5)]
write.table(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100, "AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt", sep = "\t",quote = F,append = F,col.names = F, row.names = F)




#----------- Get peaks status on CpG Islands  ----------------#
CGIs_mm10 <- read.table("CGIs_mm10.txt", header = F, stringsAsFactors = F)
head(CGIs_mm10)
CGIs_mm10_chr <- CGIs_mm10[,-c(1,7:12)]
head(CGIs_mm10_chr)
write.table(CGIs_mm10_chr, "CGIs_mm10_chr.txt", sep="\t", quote = F, append = F, row.names = F, col.names = F)

#Overlap with NRF1 extended peak center
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chr.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chr.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b CGIs_mm10_chr.txt > Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr.txt


#---------------------------  Correlation with Cofactors  -------------------------------#
#Get Cofactors data
cat jf1v2_Snp+chr.GRCm38.mm10.bed  /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38_ucsc_mm10.bed | sort -k1,1 -k2,2n > jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed
#Only the files are stored in way3_flash7 folder
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Sp1/Sp1_motif.mm10.bed  -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > SP1_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Sp1/Sp1_peaks.bed  -b SP1_motif_SNP_Indels_mm10.bed > SP1_peaks_motif_SNP_Indels_mm10.bed


bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Yy1/YY1_motif.mm10.bed  -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > YY1_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/way3_flash7/mm10/Matrix/NFY_NRF1_Peaks_features/Yy1/DMSO_YY1_rep1_peaks.bed  -b YY1_motif_SNP_Indels_mm10.bed > YY1_peaks_motif_SNP_Indels_mm10.bed

#NFY the files are taken from newway analysis
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NFY/rep1/PeakCall_R1_MACS1/NFY_rep1_peaks.bed ./
cp /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/NFY/rep2/PeakCall_R1_MACS1/NFY_rep2_peaks.bed ./
  
  
bedtools intersect -a NFY_rep1_peaks.bed -b NFY_rep2_peaks.bed  | awk '{print $1"\t"$2"\t"$3"\t"$4"_"NR"\t"$5}' > NFY_intersection_rep1_2_peaks.bed
cp /media/ankitv/Archivio2/ankit/motif_analysis/marcos/motif_trovati/mouse/motif[CCAAT].mm10.bed NFY_motif.mm10.sort.bed
bedtools intersect -wa -wb -a NFY_motif.mm10.sort.bed -b jf1v2_SNP_Grcm38.mm10_Indel.Grcm38.bed > NFY_motif_SNP_Indels_mm10.bed
bedtools intersect -wa -wb -a NFY_intersection_rep1_2_peaks.bed -b NFY_motif_SNP_Indels_mm10.bed > NFY_peaks_motif_SNP_Indels_mm10.bed


#Overlap NRF1 peaks with feature

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NRF1_B6_specific_peaks_SP1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NRF1_JF1_specific_peaks_SP1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b SP1_peaks_motif_SNP_Indels_mm10.bed > NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels.bed

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NRF1_B6_specific_peaks_YY1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NRF1_JF1_specific_peaks_YY1_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b YY1_peaks_motif_SNP_Indels_mm10.bed > NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels.bed

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt -b NFY_peaks_motif_SNP_Indels_mm10.bed > NRF1_B6_specific_peaks_NFY_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt -b NFY_peaks_motif_SNP_Indels_mm10.bed > NRF1_JF1_specific_peaks_NFY_motif_SNP_Indels.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt -b NFY_peaks_motif_SNP_Indels_mm10.bed > NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indels.bed

#------------------------ Public datasets --------------------------------#
#All files from NRF1A, NRF1B, NRF1C were combined to get merged files and lifovered to mm10
#/media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NRF1s/lift_mm10/NRF1_mm10.bed

#Overlap NRF1 
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt -b /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NRF1/lift_mm10/NRF1_mm10.bed > NRF1_B6_specific_peaks_pub_NRF1.bed

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt -b /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NRF1/lift_mm10/NRF1_mm10.bed > NRF1_JF1_specific_peaks_pub_NRF1.bed

bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b  /media/ankitv/Archivio1/download_geo_data/transcription_factors_project/mouse/chip-seq/NRF1/lift_mm10/NRF1_mm10.bed > NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1.bed


#---------------------  Countings --------------------#
#Intersections count
awk '{print $1"_"$2"_"$3"_"$4}' GCGC_jf1v2SNPmm10.bed | sort -k1,1 -u | wc -l
awk '{print $1"_"$2"_"$3"_"$4}' GCGC_jf1v2SNPmm10.jf1.bed | sort -k1,1 -u | wc -l

#NRF1_B6
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt_qual199.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt | wc -l
sort -k4,4 -u NRF1_B6_specific_peaks_GCGC.bed | wc -l
sort -k4,4 -u NRF1_B6_specific_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_B6_specific_peaks_GCGC_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_B6_specific_peaks_GCGC_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NRF1_B6_specific_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NRF1_B6_specific_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NRF1_B6_specific_peaks.bed | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_JB1InputBi.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_InputBiallelicfilt.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_midspan100.txt | wc -l

#NRF1_JF1
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt_qual199.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt | wc -l
sort -k4,4 -u NRF1_JF1_specific_peaks_GCGC.bed | wc -l
sort -k4,4 -u NRF1_JF1_specific_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_JF1_specific_peaks_GCGC_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_JF1_specific_peaks_GCGC_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NRF1_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NRF1_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NRF1_JF1_specific_peaks.bed | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_JB1InputBi.txt  | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt.txt  | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi.txt  | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt  | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_midspan100.txt | wc -l

#--------------------------------    Quality filter ------------------------------#
awk '{if($5>=199) print $0}' AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt > AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.qual199.txt
awk '{if($5>=199) print $0}' AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt > AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.qual199.txt

#-----------------------------------------  Prepare a combined Table  -------------------------------------------#
NRF1_intersection_rep1_2_peaks <- read.table("NRF1_intersection_rep1_2_peaks.bed", header = F, stringsAsFactors = F)
head(NRF1_intersection_rep1_2_peaks)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6,NRF1_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
tail(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- data.frame(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr[,c(21:23,1,24)])
colnames(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr) <- c("V1","V2","V3","V4","V5")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- data.frame(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr[order(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr$V4),]
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
head(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)
tail(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr)

AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq <- data.frame(unique(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6$V4))
colnames(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq) <- c("V4")
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq)

NRF1_B6_specific_peaks_GCGC <- read.table("NRF1_B6_specific_peaks_GCGC.bed")
NRF1_B6_specific_peaks_GCGC_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_GCGC$V4),unique(NRF1_B6_specific_peaks_GCGC$V4))
colnames(NRF1_B6_specific_peaks_GCGC_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_GCGC_uniq <- data.frame(NRF1_B6_specific_peaks_GCGC_uniq)
NRF1_B6_specific_peaks_GCGC_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_GCGC_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_GCGC_uniqB6 <- NRF1_B6_specific_peaks_GCGC_uniqB6[order(NRF1_B6_specific_peaks_GCGC_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_GCGC_uniqB6)


NRF1_B6_specific_peaks_jf1v2 <- read.table("NRF1_B6_specific_peaks_jf1v2.bed")
NRF1_B6_specific_peaks_jf1v2_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_jf1v2$V4),unique(NRF1_B6_specific_peaks_jf1v2$V4))
colnames(NRF1_B6_specific_peaks_jf1v2_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_jf1v2_uniq <- data.frame(NRF1_B6_specific_peaks_jf1v2_uniq)
NRF1_B6_specific_peaks_jf1v2_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_jf1v2_uniqB6 <- NRF1_B6_specific_peaks_jf1v2_uniqB6[order(NRF1_B6_specific_peaks_jf1v2_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_jf1v2_uniqB6)


NRF1_B6_specific_peaks_GCGC_jf1v2 <- read.table("NRF1_B6_specific_peaks_GCGC_jf1v2.bed")
NRF1_B6_specific_peaks_GCGC_jf1v2_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_GCGC_jf1v2$V4),unique(NRF1_B6_specific_peaks_GCGC_jf1v2$V4))
colnames(NRF1_B6_specific_peaks_GCGC_jf1v2_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_GCGC_jf1v2_uniq <- data.frame(NRF1_B6_specific_peaks_GCGC_jf1v2_uniq)
NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_GCGC_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6 <- NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6[order(NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6)

NRF1_B6_specific_peaks_GCGC_then_jf1v2 <- read.table("NRF1_B6_specific_peaks_GCGC_then_jf1v2.bed")
NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_GCGC_then_jf1v2$V4),unique(NRF1_B6_specific_peaks_GCGC_then_jf1v2$V4))
colnames(NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniq <- data.frame(NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniq)
NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6 <- NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6[order(NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6)#No need to include


Nomotif_NRF1_B6_specific_peaks <- read.table("Nomotif_NRF1_B6_specific_peaks.bed")
Nomotif_NRF1_B6_specific_peaks_uniq <- cbind.data.frame(unique(Nomotif_NRF1_B6_specific_peaks$V4),unique(Nomotif_NRF1_B6_specific_peaks$V4))
colnames(Nomotif_NRF1_B6_specific_peaks_uniq) <- c("V1","V4")
Nomotif_NRF1_B6_specific_peaks_uniq <- data.frame(Nomotif_NRF1_B6_specific_peaks_uniq)
Nomotif_NRF1_B6_specific_peaks_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,Nomotif_NRF1_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
Nomotif_NRF1_B6_specific_peaks_uniqB6 <- Nomotif_NRF1_B6_specific_peaks_uniqB6[order(Nomotif_NRF1_B6_specific_peaks_uniqB6$V4),]
dim(Nomotif_NRF1_B6_specific_peaks_uniqB6)


NomotifSNP_NRF1_B6_specific_peaks <- read.table("NomotifSNP_NRF1_B6_specific_peaks.bed")
NomotifSNP_NRF1_B6_specific_peaks_uniq <- cbind.data.frame(unique(NomotifSNP_NRF1_B6_specific_peaks$V4),unique(NomotifSNP_NRF1_B6_specific_peaks$V4))
colnames(NomotifSNP_NRF1_B6_specific_peaks_uniq) <- c("V1","V4")
NomotifSNP_NRF1_B6_specific_peaks_uniq <- data.frame(NomotifSNP_NRF1_B6_specific_peaks_uniq)
NomotifSNP_NRF1_B6_specific_peaks_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NomotifSNP_NRF1_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
NomotifSNP_NRF1_B6_specific_peaks_uniqB6 <- NomotifSNP_NRF1_B6_specific_peaks_uniqB6[order(NomotifSNP_NRF1_B6_specific_peaks_uniqB6$V4),]
dim(NomotifSNP_NRF1_B6_specific_peaks_uniqB6)

NoSNP_NRF1_B6_specific_peaks <- read.table("NoSNP_NRF1_B6_specific_peaks.bed")
NoSNP_NRF1_B6_specific_peaks_uniq <- cbind.data.frame(unique(NoSNP_NRF1_B6_specific_peaks$V4),unique(NoSNP_NRF1_B6_specific_peaks$V4))
colnames(NoSNP_NRF1_B6_specific_peaks_uniq) <- c("V1","V4")
NoSNP_NRF1_B6_specific_peaks_uniq <- data.frame(NoSNP_NRF1_B6_specific_peaks_uniq)
NoSNP_NRF1_B6_specific_peaks_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NoSNP_NRF1_B6_specific_peaks_uniq, by="V4", all=T,sort = F)
NoSNP_NRF1_B6_specific_peaks_uniqB6 <- NoSNP_NRF1_B6_specific_peaks_uniqB6[order(NoSNP_NRF1_B6_specific_peaks_uniqB6$V4),]
dim(NoSNP_NRF1_B6_specific_peaks_uniqB6)


AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndel <- read.table("AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndel.txt")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndel$V4),unique(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndel$V4))
colnames(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu) <- c("V1","V4")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu <- data.frame(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndelu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq[order(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq)

NRF1_B6_specific_peaks_GCGC_jf1vindel <- read.table("NRF1_B6_specific_peaks_GCGC_jf1vindel.bed")
NRF1_B6_specific_peaks_GCGC_jf1vindel_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_GCGC_jf1vindel$V4),unique(NRF1_B6_specific_peaks_GCGC_jf1vindel$V4))
colnames(NRF1_B6_specific_peaks_GCGC_jf1vindel_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_GCGC_jf1vindel_uniq <- data.frame(NRF1_B6_specific_peaks_GCGC_jf1vindel_uniq)
NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_GCGC_jf1vindel_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6 <- NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6[order(NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6)

AlspNRF1ddsNormcounts_B6_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_B6_JB1InputBi.txt")
AlspNRF1ddsNormcounts_B6_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_B6_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_B6_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_B6_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_B6_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_B6_JB1InputBiu)
AlspNRF1ddsNormcounts_B6_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_JB1InputBiuniq <- AlspNRF1ddsNormcounts_B6_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_B6_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_JB1InputBiuniq)

AlspNRF1ddsNormcounts_B6_InputBiallelicfilt <- read.table("AlspNRF1ddsNormcounts_B6_InputBiallelicfilt.txt")
AlspNRF1ddsNormcounts_B6_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt$V4),unique(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt$V4))
colnames(AlspNRF1ddsNormcounts_B6_InputBiallelicfiltu) <- c("V1","V4")
AlspNRF1ddsNormcounts_B6_InputBiallelicfiltu <- data.frame(AlspNRF1ddsNormcounts_B6_InputBiallelicfiltu)
AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq <- AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq[order(AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq)

AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi.txt")
AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiu)
AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq <- AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq)

AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR <- read.table("AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR.txt")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR$V4),unique(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICR$V4))
colnames(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu) <- c("V1","V4")
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu <- data.frame(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq[order(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq)

Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chr <- read.table("Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chr.txt")
Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chr$V4),unique(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chru <- data.frame(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chru)
Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq <- Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq[order(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq)

NRF1_B6_specific_peaks_SP1_motif_SNP_Indels <- read.table("NRF1_B6_specific_peaks_SP1_motif_SNP_Indels.bed")
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_B6_specific_peaks_SP1_motif_SNP_Indels$V4),unique(NRF1_B6_specific_peaks_SP1_motif_SNP_Indels$V4))
colnames(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsu <- data.frame(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsu)
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq[order(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq)
#Its a blank file but needed, So 
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V1 <- NA
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq <- NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq[order(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq)


NRF1_B6_specific_peaks_YY1_motif_SNP_Indels <- read.table("NRF1_B6_specific_peaks_YY1_motif_SNP_Indels.bed")
NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_B6_specific_peaks_YY1_motif_SNP_Indels$V4),unique(NRF1_B6_specific_peaks_YY1_motif_SNP_Indels$V4))
colnames(NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsu <- data.frame(NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsu)
NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq <- NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq[order(NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq)



NRF1_B6_specific_peaks_NFY_motif_SNP_Indels <- read.table("NRF1_B6_specific_peaks_NFY_motif_SNP_Indels.bed")
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_B6_specific_peaks_NFY_motif_SNP_Indels$V4),unique(NRF1_B6_specific_peaks_NFY_motif_SNP_Indels$V4))
colnames(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsu <- data.frame(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsu)
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq <- NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq[order(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq$V4),]
dim(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq)
#Its a blank file but needed, So 
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq <- AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq$V1 <- NA
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq <- NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq[order(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq$V4),]
dim(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq)



#Pub NRF1_B6
NRF1_B6_specific_peaks_pub_NRF1 <- read.table("NRF1_B6_specific_peaks_pub_NRF1.bed")
NRF1_B6_specific_peaks_pub_NRF1_uniq <- cbind.data.frame(unique(NRF1_B6_specific_peaks_pub_NRF1$V4),unique(NRF1_B6_specific_peaks_pub_NRF1$V4))
colnames(NRF1_B6_specific_peaks_pub_NRF1_uniq) <- c("V1","V4")
NRF1_B6_specific_peaks_pub_NRF1_uniq <- data.frame(NRF1_B6_specific_peaks_pub_NRF1_uniq)
NRF1_B6_specific_peaks_pub_NRF1_uniqB6 <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq,NRF1_B6_specific_peaks_pub_NRF1_uniq, by="V4", all=T,sort = F)
NRF1_B6_specific_peaks_pub_NRF1_uniqB6 <- NRF1_B6_specific_peaks_pub_NRF1_uniqB6[order(NRF1_B6_specific_peaks_pub_NRF1_uniqB6$V4),]
dim(NRF1_B6_specific_peaks_pub_NRF1_uniqB6)

#Position Het SNPs
AlspNRF1ddsNormcounts_B6_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_B6_JB1InputBi.txt")
AlspNRF1ddsNormcounts_B6_JB1InputBiSNP <- AlspNRF1ddsNormcounts_B6_JB1InputBi[,c(4,21:27)]
head(AlspNRF1ddsNormcounts_B6_JB1InputBiSNP)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNP <- data.frame(AlspNRF1ddsNormcounts_B6_JB1InputBiSNP)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq <- merge(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq, AlspNRF1ddsNormcounts_B6_JB1InputBiSNP, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq["SNP"] <- paste0(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V21,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V22,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V23,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V24,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V25,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V26,"_",
                                                           AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq <- data.frame(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)


#Assign Score 0,1
NRF1_B6_specific_peaks_GCGC_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_GCGC_uniqB6$V1),0,1)
NRF1_B6_specific_peaks_jf1v2_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_jf1v2_uniqB6$V1),0,1)
NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6$V1),0,1)
NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_GCGC_then_jf1v2_uniqB6$V1),0,1)
Nomotif_NRF1_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(Nomotif_NRF1_B6_specific_peaks_uniqB6$V1),0,1)
NomotifSNP_NRF1_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(NomotifSNP_NRF1_B6_specific_peaks_uniqB6$V1),0,1)
NoSNP_NRF1_B6_specific_peaks_uniqB6["score"] <- ifelse(is.na(NoSNP_NRF1_B6_specific_peaks_uniqB6$V1),0,1)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq$V1),0,1)
NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6$V1),0,1)
AlspNRF1ddsNormcounts_B6_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_B6_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq$V1),0,1)
AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq$V1),0,1)
Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq$V1),0,1)
NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq$V1),0,1)
NRF1_B6_specific_peaks_pub_NRF1_uniqB6["score"] <- ifelse(is.na(NRF1_B6_specific_peaks_pub_NRF1_uniqB6$V1),0,1)

matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features <- cbind.data.frame(AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_uniq_chr,
                                                                            NRF1_B6_specific_peaks_GCGC_uniqB6,
                                                                            NRF1_B6_specific_peaks_jf1v2_uniqB6,
                                                                            NRF1_B6_specific_peaks_GCGC_jf1v2_uniqB6,
                                                                            AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_jf1vIndeluniq,
                                                                            NRF1_B6_specific_peaks_GCGC_jf1vindel_uniqB6,
                                                                            AlspNRF1ddsNormcounts_B6_JB1InputBiuniq,
                                                                            AlspNRF1ddsNormcounts_B6_InputBiallelicfiltuniq,
                                                                            AlspNRF1ddsNormcounts_B6_InputBiallelicfilt_JB1InputBiuniq,
                                                                            AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_mm10_ICRuniq,
                                                                            Alsp_NRF1_B6_specific_peaks_CGIs_mm10_chruniq,
                                                                            NRF1_B6_specific_peaks_SP1_motif_SNP_Indelsuniq,
                                                                            NRF1_B6_specific_peaks_YY1_motif_SNP_Indelsuniq,
                                                                            NRF1_B6_specific_peaks_NFY_motif_SNP_Indelsuniq,
                                                                            NRF1_B6_specific_peaks_pub_NRF1_uniqB6,
                                                                            AlspNRF1ddsNormcounts_B6_JB1InputBiSNPuniq)


colnames(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features) <- c("Chr","Start","End","Peaks","Score","IDX.1","IDY.1","GCGC","IDX.2","IDY.2","SNPs","IDX.3","IDY.3","GCGC_SNPs","IDX.4","IDY.4","Indels","IDX.5","IDY.5","GCGC_Indels","IDX.6","IDY.6","Input_Het_SNPs","IDX.7","IDY.7","Input_Bi_Regions","IDX.8","IDY.8","Input_Bi_Het_SNPs","IDX.9","IDY.9","ICRs","IDX.10","IDY.10","CGIs","IDX.11","IDY.11","SP1","IDX.12","IDY.12","YY1","IDX.13","IDY.13","NFY","IDX.14","IDY.14","pub_NRF1","IDX.15","PosHetSNPs")

dim(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features)
head(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features)
tail(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features)
matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab  <- data.frame(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features[,c(1:4,seq(5,47,3),49)])
matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab["Support_Index"] <- rowMeans(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab[,6:(length(colnames(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab))-1)])
matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab <- matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab[order(-matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab)

write.table(matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab, "matrix_AlspNRF1ddsNormcounts_B6_monofilt_uniqB6_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)


#-------- NRF1_JF1
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1,NRF1_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
tail(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- data.frame(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr[,c(21:23,1,24)])
colnames(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr) <- c("V1","V2","V3","V4","V5")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- data.frame(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr[order(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr$V4),]
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
head(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)
tail(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr)

AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq <- data.frame(unique(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1$V4))
colnames(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq) <- c("V4")
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq)

NRF1_JF1_specific_peaks_GCGC <- read.table("NRF1_JF1_specific_peaks_GCGC.bed")
NRF1_JF1_specific_peaks_GCGC_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_GCGC$V4),unique(NRF1_JF1_specific_peaks_GCGC$V4))
colnames(NRF1_JF1_specific_peaks_GCGC_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_GCGC_uniq <- data.frame(NRF1_JF1_specific_peaks_GCGC_uniq)
NRF1_JF1_specific_peaks_GCGC_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_GCGC_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_GCGC_uniqJF1 <- NRF1_JF1_specific_peaks_GCGC_uniqJF1[order(NRF1_JF1_specific_peaks_GCGC_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_GCGC_uniqJF1)


NRF1_JF1_specific_peaks_jf1v2 <- read.table("NRF1_JF1_specific_peaks_jf1v2.bed")
NRF1_JF1_specific_peaks_jf1v2_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_jf1v2$V4),unique(NRF1_JF1_specific_peaks_jf1v2$V4))
colnames(NRF1_JF1_specific_peaks_jf1v2_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_jf1v2_uniq <- data.frame(NRF1_JF1_specific_peaks_jf1v2_uniq)
NRF1_JF1_specific_peaks_jf1v2_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_jf1v2_uniqJF1 <- NRF1_JF1_specific_peaks_jf1v2_uniqJF1[order(NRF1_JF1_specific_peaks_jf1v2_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_jf1v2_uniqJF1)


NRF1_JF1_specific_peaks_GCGC_jf1v2 <- read.table("NRF1_JF1_specific_peaks_GCGC_jf1v2.bed")
NRF1_JF1_specific_peaks_GCGC_jf1v2_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_GCGC_jf1v2$V4),unique(NRF1_JF1_specific_peaks_GCGC_jf1v2$V4))
colnames(NRF1_JF1_specific_peaks_GCGC_jf1v2_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_GCGC_jf1v2_uniq <- data.frame(NRF1_JF1_specific_peaks_GCGC_jf1v2_uniq)
NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_GCGC_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1 <- NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1[order(NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1)

NRF1_JF1_specific_peaks_GCGC_then_jf1v2 <- read.table("NRF1_JF1_specific_peaks_GCGC_then_jf1v2.bed")
NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_GCGC_then_jf1v2$V4),unique(NRF1_JF1_specific_peaks_GCGC_then_jf1v2$V4))
colnames(NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniq <- data.frame(NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniq)
NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1 <- NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1[order(NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1)#No need to include


Nomotif_NRF1_JF1_specific_peaks <- read.table("Nomotif_NRF1_JF1_specific_peaks.bed")
Nomotif_NRF1_JF1_specific_peaks_uniq <- cbind.data.frame(unique(Nomotif_NRF1_JF1_specific_peaks$V4),unique(Nomotif_NRF1_JF1_specific_peaks$V4))
colnames(Nomotif_NRF1_JF1_specific_peaks_uniq) <- c("V1","V4")
Nomotif_NRF1_JF1_specific_peaks_uniq <- data.frame(Nomotif_NRF1_JF1_specific_peaks_uniq)
Nomotif_NRF1_JF1_specific_peaks_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,Nomotif_NRF1_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
Nomotif_NRF1_JF1_specific_peaks_uniqJF1 <- Nomotif_NRF1_JF1_specific_peaks_uniqJF1[order(Nomotif_NRF1_JF1_specific_peaks_uniqJF1$V4),]
dim(Nomotif_NRF1_JF1_specific_peaks_uniqJF1)


NomotifSNP_NRF1_JF1_specific_peaks <- read.table("NomotifSNP_NRF1_JF1_specific_peaks.bed")
NomotifSNP_NRF1_JF1_specific_peaks_uniq <- cbind.data.frame(unique(NomotifSNP_NRF1_JF1_specific_peaks$V4),unique(NomotifSNP_NRF1_JF1_specific_peaks$V4))
colnames(NomotifSNP_NRF1_JF1_specific_peaks_uniq) <- c("V1","V4")
NomotifSNP_NRF1_JF1_specific_peaks_uniq <- data.frame(NomotifSNP_NRF1_JF1_specific_peaks_uniq)
NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NomotifSNP_NRF1_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1 <- NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1[order(NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1$V4),]
dim(NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1)

NoSNP_NRF1_JF1_specific_peaks <- read.table("NoSNP_NRF1_JF1_specific_peaks.bed")
NoSNP_NRF1_JF1_specific_peaks_uniq <- cbind.data.frame(unique(NoSNP_NRF1_JF1_specific_peaks$V4),unique(NoSNP_NRF1_JF1_specific_peaks$V4))
colnames(NoSNP_NRF1_JF1_specific_peaks_uniq) <- c("V1","V4")
NoSNP_NRF1_JF1_specific_peaks_uniq <- data.frame(NoSNP_NRF1_JF1_specific_peaks_uniq)
NoSNP_NRF1_JF1_specific_peaks_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NoSNP_NRF1_JF1_specific_peaks_uniq, by="V4", all=T,sort = F)
NoSNP_NRF1_JF1_specific_peaks_uniqJF1 <- NoSNP_NRF1_JF1_specific_peaks_uniqJF1[order(NoSNP_NRF1_JF1_specific_peaks_uniqJF1$V4),]
dim(NoSNP_NRF1_JF1_specific_peaks_uniqJF1)

AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel <- read.table("AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel.txt")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel$V4),unique(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndel$V4))
colnames(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu) <- c("V1","V4")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu <- data.frame(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndelu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq[order(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq)

NRF1_JF1_specific_peaks_GCGC_jf1vindel <- read.table("NRF1_JF1_specific_peaks_GCGC_jf1vindel.bed")
NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_GCGC_jf1vindel$V4),unique(NRF1_JF1_specific_peaks_GCGC_jf1vindel$V4))
colnames(NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniq <- data.frame(NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniq)
NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1 <- NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1[order(NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1)

AlspNRF1ddsNormcounts_JF1_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_JF1_JB1InputBi.txt")
AlspNRF1ddsNormcounts_JF1_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_JF1_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_JF1_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_JF1_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_JF1_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_JF1_JB1InputBiu)
AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq <- AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq)

AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt <- read.table("AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt.txt")
AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt$V4),unique(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt$V4))
colnames(AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltu) <- c("V1","V4")
AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltu <- data.frame(AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltu)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq <- AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq[order(AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq)

AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi.txt")
AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiu)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq <- AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq)

AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR <- read.table("AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR.txt")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR$V4),unique(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICR$V4))
colnames(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu) <- c("V1","V4")
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu <- data.frame(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq[order(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq)
#Its a blank file but needed, So 
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V1 <- NA
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq[order(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq)

Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chr <- read.table("Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chr.txt")
Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chr$V4),unique(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chru <- data.frame(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chru)
Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq <- Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq[order(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq)

NRF1_JF1_specific_peaks_SP1_motif_SNP_Indels <- read.table("NRF1_JF1_specific_peaks_SP1_motif_SNP_Indels.bed")
NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indels$V4),unique(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indels$V4))
colnames(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsu <- data.frame(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsu)
NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq <- NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq[order(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq)

NRF1_JF1_specific_peaks_YY1_motif_SNP_Indels <- read.table("NRF1_JF1_specific_peaks_YY1_motif_SNP_Indels.bed")
NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indels$V4),unique(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indels$V4))
colnames(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsu <- data.frame(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsu)
NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq <- NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq[order(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq)

NRF1_JF1_specific_peaks_NFY_motif_SNP_Indels <- read.table("NRF1_JF1_specific_peaks_NFY_motif_SNP_Indels.bed")
NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indels$V4),unique(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indels$V4))
colnames(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsu <- data.frame(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsu)
NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq <- NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq[order(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq$V4),]
dim(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq)


#Pub NRF1_JF1
NRF1_JF1_specific_peaks_pub_NRF1 <- read.table("NRF1_JF1_specific_peaks_pub_NRF1.bed")
NRF1_JF1_specific_peaks_pub_NRF1_uniq <- cbind.data.frame(unique(NRF1_JF1_specific_peaks_pub_NRF1$V4),unique(NRF1_JF1_specific_peaks_pub_NRF1$V4))
colnames(NRF1_JF1_specific_peaks_pub_NRF1_uniq) <- c("V1","V4")
NRF1_JF1_specific_peaks_pub_NRF1_uniq <- data.frame(NRF1_JF1_specific_peaks_pub_NRF1_uniq)
NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1 <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq,NRF1_JF1_specific_peaks_pub_NRF1_uniq, by="V4", all=T,sort = F)
NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1 <- NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1[order(NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1$V4),]
dim(NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1)


#Position Het SNPs
AlspNRF1ddsNormcounts_JF1_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_JF1_JB1InputBi.txt")
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNP <- AlspNRF1ddsNormcounts_JF1_JB1InputBi[,c(4,21:27)]
head(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNP)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNP <- data.frame(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNP)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq <- merge(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq, AlspNRF1ddsNormcounts_JF1_JB1InputBiSNP, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq["SNP"] <- paste0(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V21,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V22,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V23,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V24,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V25,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V26,"_",
                                                            AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq <- data.frame(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)

NRF1_JF1_specific_peaks_GCGC_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_GCGC_uniqJF1$V1),0,1)
NRF1_JF1_specific_peaks_jf1v2_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_jf1v2_uniqJF1$V1),0,1)
NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1$V1),0,1)
NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_GCGC_then_jf1v2_uniqJF1$V1),0,1)
Nomotif_NRF1_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(Nomotif_NRF1_JF1_specific_peaks_uniqJF1$V1),0,1)
NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(NomotifSNP_NRF1_JF1_specific_peaks_uniqJF1$V1),0,1)
NoSNP_NRF1_JF1_specific_peaks_uniqJF1["score"] <- ifelse(is.na(NoSNP_NRF1_JF1_specific_peaks_uniqJF1$V1),0,1)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq$V1),0,1)
NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1$V1),0,1)
AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq$V1),0,1)
AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq$V1),0,1)
Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq$V1),0,1)
NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq$V1),0,1)
NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1["score"] <- ifelse(is.na(NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1$V1),0,1)


matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features <- cbind.data.frame(AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_uniq_chr,
                                                                              NRF1_JF1_specific_peaks_GCGC_uniqJF1,
                                                                              NRF1_JF1_specific_peaks_jf1v2_uniqJF1,
                                                                              NRF1_JF1_specific_peaks_GCGC_jf1v2_uniqJF1,
                                                                              AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_jf1vIndeluniq,
                                                                              NRF1_JF1_specific_peaks_GCGC_jf1vindel_uniqJF1,
                                                                              AlspNRF1ddsNormcounts_JF1_JB1InputBiuniq,
                                                                              AlspNRF1ddsNormcounts_JF1_InputBiallelicfiltuniq,
                                                                              AlspNRF1ddsNormcounts_JF1_InputBiallelicfilt_JB1InputBiuniq,
                                                                              AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_mm10_ICRuniq,
                                                                              Alsp_NRF1_JF1_specific_peaks_CGIs_mm10_chruniq,
                                                                              NRF1_JF1_specific_peaks_SP1_motif_SNP_Indelsuniq,
                                                                              NRF1_JF1_specific_peaks_YY1_motif_SNP_Indelsuniq,
                                                                              NRF1_JF1_specific_peaks_NFY_motif_SNP_Indelsuniq,
                                                                              NRF1_JF1_specific_peaks_pub_NRF1_uniqJF1,
                                                                              AlspNRF1ddsNormcounts_JF1_JB1InputBiSNPuniq)





colnames(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features) <- c("Chr","Start","End","Peaks","Score","IDX.1","IDY.1","GCGC","IDX.2","IDY.2","SNPs","IDX.3","IDY.3","GCGC_SNPs","IDX.4","IDY.4","Indels","IDX.5","IDY.5","GCGC_Indels","IDX.6","IDY.6","Input_Het_SNPs","IDX.7","IDY.7","Input_Bi_Regions","IDX.8","IDY.8","Input_Bi_Het_SNPs","IDX.9","IDY.9","ICRs","IDX.10","IDY.10","CGIs","IDX.11","IDY.11","SP1","IDX.12","IDY.12","YY1","IDX.13","IDY.13","NFY","IDX.14","IDY.14","pub_NRF1","IDX.15","PosHetSNPs")

dim(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features)
head(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features)
tail(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features)
matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab  <- data.frame(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features[,c(1:4,seq(5,47,3),49)])
matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab["Support_Index"] <- rowMeans(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab[,6:(length(colnames(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab))-1)])
matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab <- matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab[order(-matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab)

write.table(matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab, "matrix_AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)



#----------------------------- Biallelic peaks intersection ------------------
#AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt

#NRF1
#Motif and SNP overlap peaks
#GCGC
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_motif.mm10.sort.bed > NRF1_B_biallelic_peaks_GCGC.bed
#SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_B_biallelic_peaks_jf1v2.bed

#Note GCGC_jf1v2SNPmm10.bed is file intersected with SNP (jf1v2mm10)
#GCGC+SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1v2SNPmm10.bed > NRF1_B_biallelic_peaks_GCGC_jf1v2.bed

#GCGC, SNP
bedtools intersect -wa -wb -a NRF1_B_biallelic_peaks_GCGC.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_B_biallelic_peaks_GCGC_then_jf1v2.bed

#NO NRF1
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_motif.mm10.sort.bed -v > Nomotif_NRF1_B_biallelic_peaks.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1v2SNPmm10.bed -v > NomotifSNP_NRF1_B_biallelic_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NRF1_B_biallelic_peaks.bed

#Motif and SNP overlap peaks JF1
#GCGC
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_motif.mm10.jf1.sort.bed > NRF1_J_biallelic_peaks_GCGC.bed
#SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_J_biallelic_peaks_jf1v2.bed
#Note GCGC_jf1v2SNPmm10.jf1.bed is file intersected with SNP (jf1v2mm10.jf1)
#GCGC+SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1v2SNPmm10.jf1.bed > NRF1_J_biallelic_peaks_GCGC_jf1v2.bed
#GCGC, SNP
bedtools intersect -wa -wb -a NRF1_J_biallelic_peaks_GCGC.bed -b jf1v2_Snp+chr.GRCm38.mm10.bed > NRF1_J_biallelic_peaks_GCGC_then_jf1v2.bed

#NO NRF1
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_motif.mm10.jf1.sort.bed -v > Nomotif_NRF1_J_biallelic_peaks.bed
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1v2SNPmm10.jf1.bed -v > NomotifSNP_NRF1_J_biallelic_peaks.bed

#No SNP
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b jf1v2_Snp+chr.GRCm38.mm10.bed -v > NoSNP_NRF1_J_biallelic_peaks.bed


bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b JB1resultInputBi_chr.txt > AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b AlspInputddsNormcounts_Biallelicfilt.txt > AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b AlspInputddsNormcounts_Biallelicfilt_JB1resultInputBi.txt > AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi.txt
bedtools intersect -wa -wb -a  AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt  -b AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt > AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt

#Indel
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b /media/ankitv/Archivio1/MEA/test/mm10/bcftools/JF1_Indel_Grcm38.bed  > AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt

#Motif, Indels, as before
#GCGC+Indels,B6 based
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1vindelSNPmm10.bed > NRF1_B_biallelic_peaks_GCGC_jf1vindel.bed

#GCGC+Indels,JF1-based
bedtools intersect -wa -wb -a AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt -b GCGC_jf1vindelSNPmm10.jf1.bed > NRF1_J_biallelic_peaks_GCGC_jf1vindel.bed


#Countings
#Common NRF1 Biallelic B and J
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_qual199.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt | wc -l
sort -k4,4 -u AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_midspan100.txt | wc -l

#NRF1_B
sort -k4,4 -u NRF1_B_biallelic_peaks_GCGC_jf1vindel.bed | wc -l
sort -k4,4 -u NRF1_B_biallelic_peaks_GCGC.bed | wc -l
sort -k4,4 -u NRF1_B_biallelic_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_B_biallelic_peaks_GCGC_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_B_biallelic_peaks_GCGC_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NRF1_B_biallelic_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NRF1_B_biallelic_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NRF1_B_biallelic_peaks.bed | wc -l

#NRF1_J
sort -k4,4 -u NRF1_J_biallelic_peaks_GCGC_jf1vindel.bed | wc -l
sort -k4,4 -u NRF1_J_biallelic_peaks_GCGC.bed | wc -l
sort -k4,4 -u NRF1_J_biallelic_peaks_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_J_biallelic_peaks_GCGC_jf1v2.bed | wc -l
sort -k4,4 -u NRF1_J_biallelic_peaks_GCGC_then_jf1v2.bed | wc -l
sort -k4,4 -u Nomotif_NRF1_J_biallelic_peaks.bed | wc -l
sort -k4,4 -u NomotifSNP_NRF1_J_biallelic_peaks.bed | wc -l
sort -k4,4 -u NoSNP_NRF1_J_biallelic_peaks.bed | wc -l


awk '{if($5>=199) print $0}' AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.txt > AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1.qual199.txt

##-- Prepare Table NRF1 Biallelic Peaks --##
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1,NRF1_intersection_rep1_2_peaks, by="V4", all.x=T,all.y=F,sort = F)
head(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
dim(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
tail(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr[,c(21:23,1,24)])
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr) <- c("V1","V2","V3","V4","V5")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr[order(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
head(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)
tail(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr)

AlspNRF1ddsNormcounts_Biallelicfilt_shared <- data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_shared) <- c("V4")
dim(AlspNRF1ddsNormcounts_Biallelicfilt_shared)


NRF1_B_biallelic_peaks_GCGC <- read.table("NRF1_B_biallelic_peaks_GCGC.bed")
NRF1_B_biallelic_peaks_GCGC_shared <- cbind.data.frame(unique(NRF1_B_biallelic_peaks_GCGC$V4),unique(NRF1_B_biallelic_peaks_GCGC$V4))
colnames(NRF1_B_biallelic_peaks_GCGC_shared) <- c("V1","V4")
NRF1_B_biallelic_peaks_GCGC_shared <- data.frame(NRF1_B_biallelic_peaks_GCGC_shared)
NRF1_B_biallelic_peaks_GCGC_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_B_biallelic_peaks_GCGC_shared, by="V4", all=T,sort = F)
NRF1_B_biallelic_peaks_GCGC_shared <- NRF1_B_biallelic_peaks_GCGC_shared[order(NRF1_B_biallelic_peaks_GCGC_shared$V4),]
dim(NRF1_B_biallelic_peaks_GCGC_shared)


NRF1_B_biallelic_peaks_jf1v2 <- read.table("NRF1_B_biallelic_peaks_jf1v2.bed")
NRF1_B_biallelic_peaks_jf1v2_shared <- cbind.data.frame(unique(NRF1_B_biallelic_peaks_jf1v2$V4),unique(NRF1_B_biallelic_peaks_jf1v2$V4))
colnames(NRF1_B_biallelic_peaks_jf1v2_shared) <- c("V1","V4")
NRF1_B_biallelic_peaks_jf1v2_shared <- data.frame(NRF1_B_biallelic_peaks_jf1v2_shared)
NRF1_B_biallelic_peaks_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_B_biallelic_peaks_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_B_biallelic_peaks_jf1v2_shared <- NRF1_B_biallelic_peaks_jf1v2_shared[order(NRF1_B_biallelic_peaks_jf1v2_shared$V4),]
dim(NRF1_B_biallelic_peaks_jf1v2_shared)


NRF1_B_biallelic_peaks_GCGC_jf1v2 <- read.table("NRF1_B_biallelic_peaks_GCGC_jf1v2.bed")
NRF1_B_biallelic_peaks_GCGC_jf1v2_shared <- cbind.data.frame(unique(NRF1_B_biallelic_peaks_GCGC_jf1v2$V4),unique(NRF1_B_biallelic_peaks_GCGC_jf1v2$V4))
colnames(NRF1_B_biallelic_peaks_GCGC_jf1v2_shared) <- c("V1","V4")
NRF1_B_biallelic_peaks_GCGC_jf1v2_shared <- data.frame(NRF1_B_biallelic_peaks_GCGC_jf1v2_shared)
NRF1_B_biallelic_peaks_GCGC_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_B_biallelic_peaks_GCGC_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_B_biallelic_peaks_GCGC_jf1v2_shared <- NRF1_B_biallelic_peaks_GCGC_jf1v2_shared[order(NRF1_B_biallelic_peaks_GCGC_jf1v2_shared$V4),]
dim(NRF1_B_biallelic_peaks_GCGC_jf1v2_shared)

NRF1_B_biallelic_peaks_GCGC_then_jf1v2 <- read.table("NRF1_B_biallelic_peaks_GCGC_then_jf1v2.bed")
NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared <- cbind.data.frame(unique(NRF1_B_biallelic_peaks_GCGC_then_jf1v2$V4),unique(NRF1_B_biallelic_peaks_GCGC_then_jf1v2$V4))
colnames(NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared) <- c("V1","V4")
NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared <- data.frame(NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared)
NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared <- NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared[order(NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared$V4),]
dim(NRF1_B_biallelic_peaks_GCGC_then_jf1v2_shared)#No need to include


Nomotif_NRF1_B_biallelic_peaks <- read.table("Nomotif_NRF1_B_biallelic_peaks.bed")
Nomotif_NRF1_B_biallelic_peaks_shared <- cbind.data.frame(unique(Nomotif_NRF1_B_biallelic_peaks$V4),unique(Nomotif_NRF1_B_biallelic_peaks$V4))
colnames(Nomotif_NRF1_B_biallelic_peaks_shared) <- c("V1","V4")
Nomotif_NRF1_B_biallelic_peaks_shared <- data.frame(Nomotif_NRF1_B_biallelic_peaks_shared)
Nomotif_NRF1_B_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,Nomotif_NRF1_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
Nomotif_NRF1_B_biallelic_peaks_shared <- Nomotif_NRF1_B_biallelic_peaks_shared[order(Nomotif_NRF1_B_biallelic_peaks_shared$V4),]
dim(Nomotif_NRF1_B_biallelic_peaks_shared)


NomotifSNP_NRF1_B_biallelic_peaks <- read.table("NomotifSNP_NRF1_B_biallelic_peaks.bed")
NomotifSNP_NRF1_B_biallelic_peaks_shared <- cbind.data.frame(unique(NomotifSNP_NRF1_B_biallelic_peaks$V4),unique(NomotifSNP_NRF1_B_biallelic_peaks$V4))
colnames(NomotifSNP_NRF1_B_biallelic_peaks_shared) <- c("V1","V4")
NomotifSNP_NRF1_B_biallelic_peaks_shared <- data.frame(NomotifSNP_NRF1_B_biallelic_peaks_shared)
NomotifSNP_NRF1_B_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NomotifSNP_NRF1_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
NomotifSNP_NRF1_B_biallelic_peaks_shared <- NomotifSNP_NRF1_B_biallelic_peaks_shared[order(NomotifSNP_NRF1_B_biallelic_peaks_shared$V4),]
dim(NomotifSNP_NRF1_B_biallelic_peaks_shared)

NoSNP_NRF1_B_biallelic_peaks <- read.table("NoSNP_NRF1_B_biallelic_peaks.bed")
NoSNP_NRF1_B_biallelic_peaks_shared <- cbind.data.frame(unique(NoSNP_NRF1_B_biallelic_peaks$V4),unique(NoSNP_NRF1_B_biallelic_peaks$V4))
colnames(NoSNP_NRF1_B_biallelic_peaks_shared) <- c("V1","V4")
NoSNP_NRF1_B_biallelic_peaks_shared <- data.frame(NoSNP_NRF1_B_biallelic_peaks_shared)
NoSNP_NRF1_B_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NoSNP_NRF1_B_biallelic_peaks_shared, by="V4", all=T,sort = F)
NoSNP_NRF1_B_biallelic_peaks_shared <- NoSNP_NRF1_B_biallelic_peaks_shared[order(NoSNP_NRF1_B_biallelic_peaks_shared$V4),]
dim(NoSNP_NRF1_B_biallelic_peaks_shared)

NRF1_B_biallelic_peaks_GCGC_jf1vindel <- read.table("NRF1_B_biallelic_peaks_GCGC_jf1vindel.bed")
NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniq <- cbind.data.frame(unique(NRF1_B_biallelic_peaks_GCGC_jf1vindel$V4),unique(NRF1_B_biallelic_peaks_GCGC_jf1vindel$V4))
colnames(NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniq) <- c("V1","V4")
NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniq <- data.frame(NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniq)
NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6 <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared,NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniq, by="V4", all=T,sort = F)
NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6 <- NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6[order(NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6$V4),]
dim(NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6)



NRF1_J_biallelic_peaks_GCGC <- read.table("NRF1_J_biallelic_peaks_GCGC.bed")
NRF1_J_biallelic_peaks_GCGC_shared <- cbind.data.frame(unique(NRF1_J_biallelic_peaks_GCGC$V4),unique(NRF1_J_biallelic_peaks_GCGC$V4))
colnames(NRF1_J_biallelic_peaks_GCGC_shared) <- c("V1","V4")
NRF1_J_biallelic_peaks_GCGC_shared <- data.frame(NRF1_J_biallelic_peaks_GCGC_shared)
NRF1_J_biallelic_peaks_GCGC_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_J_biallelic_peaks_GCGC_shared, by="V4", all=T,sort = F)
NRF1_J_biallelic_peaks_GCGC_shared <- NRF1_J_biallelic_peaks_GCGC_shared[order(NRF1_J_biallelic_peaks_GCGC_shared$V4),]
dim(NRF1_J_biallelic_peaks_GCGC_shared)


NRF1_J_biallelic_peaks_jf1v2 <- read.table("NRF1_J_biallelic_peaks_jf1v2.bed")
NRF1_J_biallelic_peaks_jf1v2_shared <- cbind.data.frame(unique(NRF1_J_biallelic_peaks_jf1v2$V4),unique(NRF1_J_biallelic_peaks_jf1v2$V4))
colnames(NRF1_J_biallelic_peaks_jf1v2_shared) <- c("V1","V4")
NRF1_J_biallelic_peaks_jf1v2_shared <- data.frame(NRF1_J_biallelic_peaks_jf1v2_shared)
NRF1_J_biallelic_peaks_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_J_biallelic_peaks_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_J_biallelic_peaks_jf1v2_shared <- NRF1_J_biallelic_peaks_jf1v2_shared[order(NRF1_J_biallelic_peaks_jf1v2_shared$V4),]
dim(NRF1_J_biallelic_peaks_jf1v2_shared)


NRF1_J_biallelic_peaks_GCGC_jf1v2 <- read.table("NRF1_J_biallelic_peaks_GCGC_jf1v2.bed")
NRF1_J_biallelic_peaks_GCGC_jf1v2_shared <- cbind.data.frame(unique(NRF1_J_biallelic_peaks_GCGC_jf1v2$V4),unique(NRF1_J_biallelic_peaks_GCGC_jf1v2$V4))
colnames(NRF1_J_biallelic_peaks_GCGC_jf1v2_shared) <- c("V1","V4")
NRF1_J_biallelic_peaks_GCGC_jf1v2_shared <- data.frame(NRF1_J_biallelic_peaks_GCGC_jf1v2_shared)
NRF1_J_biallelic_peaks_GCGC_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_J_biallelic_peaks_GCGC_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_J_biallelic_peaks_GCGC_jf1v2_shared <- NRF1_J_biallelic_peaks_GCGC_jf1v2_shared[order(NRF1_J_biallelic_peaks_GCGC_jf1v2_shared$V4),]
dim(NRF1_J_biallelic_peaks_GCGC_jf1v2_shared)

NRF1_J_biallelic_peaks_GCGC_then_jf1v2 <- read.table("NRF1_J_biallelic_peaks_GCGC_then_jf1v2.bed")
NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared <- cbind.data.frame(unique(NRF1_J_biallelic_peaks_GCGC_then_jf1v2$V4),unique(NRF1_J_biallelic_peaks_GCGC_then_jf1v2$V4))
colnames(NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared) <- c("V1","V4")
NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared <- data.frame(NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared)
NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared, by="V4", all=T,sort = F)
NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared <- NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared[order(NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared$V4),]
dim(NRF1_J_biallelic_peaks_GCGC_then_jf1v2_shared)#No need to include


Nomotif_NRF1_J_biallelic_peaks <- read.table("Nomotif_NRF1_J_biallelic_peaks.bed")
Nomotif_NRF1_J_biallelic_peaks_shared <- cbind.data.frame(unique(Nomotif_NRF1_J_biallelic_peaks$V4),unique(Nomotif_NRF1_J_biallelic_peaks$V4))
colnames(Nomotif_NRF1_J_biallelic_peaks_shared) <- c("V1","V4")
Nomotif_NRF1_J_biallelic_peaks_shared <- data.frame(Nomotif_NRF1_J_biallelic_peaks_shared)
Nomotif_NRF1_J_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,Nomotif_NRF1_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
Nomotif_NRF1_J_biallelic_peaks_shared <- Nomotif_NRF1_J_biallelic_peaks_shared[order(Nomotif_NRF1_J_biallelic_peaks_shared$V4),]
dim(Nomotif_NRF1_J_biallelic_peaks_shared)


NomotifSNP_NRF1_J_biallelic_peaks <- read.table("NomotifSNP_NRF1_J_biallelic_peaks.bed")
NomotifSNP_NRF1_J_biallelic_peaks_shared <- cbind.data.frame(unique(NomotifSNP_NRF1_J_biallelic_peaks$V4),unique(NomotifSNP_NRF1_J_biallelic_peaks$V4))
colnames(NomotifSNP_NRF1_J_biallelic_peaks_shared) <- c("V1","V4")
NomotifSNP_NRF1_J_biallelic_peaks_shared <- data.frame(NomotifSNP_NRF1_J_biallelic_peaks_shared)
NomotifSNP_NRF1_J_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NomotifSNP_NRF1_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
NomotifSNP_NRF1_J_biallelic_peaks_shared <- NomotifSNP_NRF1_J_biallelic_peaks_shared[order(NomotifSNP_NRF1_J_biallelic_peaks_shared$V4),]
dim(NomotifSNP_NRF1_J_biallelic_peaks_shared)

NoSNP_NRF1_J_biallelic_peaks <- read.table("NoSNP_NRF1_J_biallelic_peaks.bed")
NoSNP_NRF1_J_biallelic_peaks_shared <- cbind.data.frame(unique(NoSNP_NRF1_J_biallelic_peaks$V4),unique(NoSNP_NRF1_J_biallelic_peaks$V4))
colnames(NoSNP_NRF1_J_biallelic_peaks_shared) <- c("V1","V4")
NoSNP_NRF1_J_biallelic_peaks_shared <- data.frame(NoSNP_NRF1_J_biallelic_peaks_shared)
NoSNP_NRF1_J_biallelic_peaks_shared <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared ,NoSNP_NRF1_J_biallelic_peaks_shared, by="V4", all=T,sort = F)
NoSNP_NRF1_J_biallelic_peaks_shared <- NoSNP_NRF1_J_biallelic_peaks_shared[order(NoSNP_NRF1_J_biallelic_peaks_shared$V4),]
dim(NoSNP_NRF1_J_biallelic_peaks_shared)

AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel$V4),unique(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndel$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu) <- c("V1","V4")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndelu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq)

NRF1_J_biallelic_peaks_GCGC_jf1vindel <- read.table("NRF1_J_biallelic_peaks_GCGC_jf1vindel.bed")
NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniq <- cbind.data.frame(unique(NRF1_J_biallelic_peaks_GCGC_jf1vindel$V4),unique(NRF1_J_biallelic_peaks_GCGC_jf1vindel$V4))
colnames(NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniq) <- c("V1","V4")
NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniq <- data.frame(NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniq)
NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1 <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared,NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniq, by="V4", all=T,sort = F)
NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1 <- NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1[order(NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1$V4),]
dim(NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1)

AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiu)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq)

AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt$V4),unique(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltu) <- c("V1","V4")
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltu <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltu)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq)

AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi$V4),unique(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBi$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiu) <- c("V1","V4")
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiu <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiu)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq)

AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu <- cbind.data.frame(unique(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR$V4),unique(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICR$V4))
colnames(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu) <- c("V1","V4")
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRu, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq)

#Its a blank file but needed, So 
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_shared
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V1 <- NA
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq)


Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr <- read.table("Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr.txt")
Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru <- cbind.data.frame(unique(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr$V4),unique(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chr$V4))
colnames(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru) <- c("V1","V4")
Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru <- data.frame(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru)
Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chru, by="V4", all=T,sort = F)
Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq <- Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq[order(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq$V4),]
dim(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq)


NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels <- read.table("NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels.bed")
NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels$V4),unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indels$V4))
colnames(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu <- data.frame(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu)
NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq <- NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq[order(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq)

NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels <- read.table("NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels.bed")
NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels$V4),unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indels$V4))
colnames(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu <- data.frame(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu)
NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq <- NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq[order(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq$V4),]
dim(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq)


NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indels <- read.table("NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indels.bed")
NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsu <- cbind.data.frame(unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indels$V4),unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indels$V4))
colnames(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsu) <- c("V1","V4")
NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsu <- data.frame(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsu)
NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsu, by="V4", all=T,sort = F)
NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq <- NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq[order(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq$V4),]
dim(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq)

#Pub NRF1_Biallelic
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1 <- read.table("NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1.bed")
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniq <- cbind.data.frame(unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1$V4),unique(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1$V4))
colnames(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniq) <- c("V1","V4")
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniq <- data.frame(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniq)
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared,NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniq, by="V4", all=T,sort = F)
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi <- NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi[order(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi$V4),]
dim(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi)



#Position Het SNPs
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi <- read.table("AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi.txt")
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNP <- AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBi[,c(4,21:27)]
head(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNP)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNP <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNP)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq <- merge(AlspNRF1ddsNormcounts_Biallelicfilt_shared, AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNP, by="V4", all=T,sort = F)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq["SNP"] <- paste0(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V21,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V22,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V23,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V24,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V25,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V26,"_",
                                                                      AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V27)
library(tidyverse)
library(stringr)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq %>%
  group_by(V4) %>%
  summarize(SNPs = str_c(SNP, collapse = ","))
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq <- data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq <- AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq[order(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq$V4),]
dim(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)
head(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)



#Assign Score 0,1
NRF1_B_biallelic_peaks_GCGC_shared["score"] <- ifelse(is.na(NRF1_B_biallelic_peaks_GCGC_shared$V1),0,1)
NRF1_B_biallelic_peaks_jf1v2_shared["score"] <- ifelse(is.na(NRF1_B_biallelic_peaks_jf1v2_shared$V1),0,1)
NRF1_B_biallelic_peaks_GCGC_jf1v2_shared["score"] <- ifelse(is.na(NRF1_B_biallelic_peaks_GCGC_jf1v2_shared$V1),0,1)
NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6["score"] <- ifelse(is.na(NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6$V1),0,1)
NRF1_J_biallelic_peaks_GCGC_shared["score"] <- ifelse(is.na(NRF1_J_biallelic_peaks_GCGC_shared$V1),0,1)
NRF1_J_biallelic_peaks_jf1v2_shared["score"] <- ifelse(is.na(NRF1_J_biallelic_peaks_jf1v2_shared$V1),0,1)
NRF1_J_biallelic_peaks_GCGC_jf1v2_shared["score"] <- ifelse(is.na(NRF1_J_biallelic_peaks_GCGC_jf1v2_shared$V1),0,1)
NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1["score"] <- ifelse(is.na(NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1$V1),0,1)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq$V1),0,1)
AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq$V1),0,1)
AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq$V1),0,1)
AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq["score"] <- ifelse(is.na(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq$V1),0,1)
Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq["score"] <- ifelse(is.na(Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq$V1),0,1)
NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq$V1),0,1)
NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq["score"] <- ifelse(is.na(NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq$V1),0,1)
NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi["score"] <- ifelse(is.na(NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi$V1),0,1)


matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features <- cbind.data.frame(AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_shared_chr,
                                                                       NRF1_B_biallelic_peaks_GCGC_shared,
                                                                       NRF1_B_biallelic_peaks_jf1v2_shared,
                                                                       NRF1_B_biallelic_peaks_GCGC_jf1v2_shared,
                                                                       NRF1_B_biallelic_peaks_GCGC_jf1vindel_uniqB6,
                                                                       NRF1_J_biallelic_peaks_GCGC_shared,
                                                                       NRF1_J_biallelic_peaks_jf1v2_shared,
                                                                       NRF1_J_biallelic_peaks_GCGC_jf1v2_shared,
                                                                       NRF1_J_biallelic_peaks_GCGC_jf1vindel_uniqJF1,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_jf1vIndeluniq,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiuniq,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfiltuniq,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_InputBiallelicfilt_JB1InputBiuniq,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_sharedB6JF1_mm10_ICRuniq,
                                                                       Alsp_NRF1_Biallelicfilt_sharedB6JF1_peaks_CGIs_mm10_chruniq,
                                                                       NRF1_Biallelicfilt_sharedB6JF1_peaks_SP1_motif_SNP_Indelsuniq,
                                                                       NRF1_Biallelicfilt_sharedB6JF1_peaks_YY1_motif_SNP_Indelsuniq,
                                                                       NRF1_Biallelicfilt_sharedB6JF1_peaks_NFY_motif_SNP_Indelsuniq,
                                                                       NRF1_Biallelicfilt_sharedB6JF1_peaks_pub_NRF1_uniqBi,
                                                                       AlspNRF1ddsNormcounts_Biallelicfilt_JB1InputBiSNPuniq)

colnames(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features) <- c("Chr","Start","End","PeakID","Score","IDX.1","IDY.1","B_GCGC","IDX.2","IDY.2","B_SNPs","IDX.3","IDY.3","B_GCGC_SNPs","IDX.4","IDY.4","B_GCGC_Indel","IDX.5","IDY.5","J_GCGC","IDX.6","IDY.6","J_SNPs","IDX.7","IDY.7","J_GCGC_SNPs","IDX.8","IDY.8","J_GCGC_Indel","IDX.9","IDY.9","sharedB6JF1_Indel","IDX.10","IDY.10","Input_Het_SNPs","IDX.11","IDY.11","Input_Bi_Regions","IDX.12","IDY.12","Input_Bi_Het_SNPs","IDX.13","IDY.13","ICRs","IDX.14","IDY.14","CGIs","IDX.15","IDY.15","SP1","IDX.16","IDY.16","YY1","IDX.17","IDY.17","NFY","IDX.18","IDY.18","pub_uniqBi","IDX.19","PosHetSNPs")

dim(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features)
head(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features)
tail(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features)
matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab  <- data.frame(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features[,c(1:4,seq(5,59,3),61)])
matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab["Support_Index"] <- rowMeans(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab[,6:(length(colnames(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab))-1)])
matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab <- matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab[order(-matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab$Input_Bi_Het_SNPs),]
head(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab)

write.table(matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab, "matrix_AlspNRF1ddsNormcounts_Biallelicfilt_features_tab.txt", sep="\t", quote=F, append = F, row.names = F)


#Upload reads for visulaization
sort -k1,1 -k2,2n /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_1_B6_reads.sort.bed
sort -k1,1 -k2,2n  /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed > Input_2_B6_reads.sort.bed
sort -k1,1 -k2,2n  /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_1_JF1_reads.sort.bed
sort -k1,1 -k2,2n  /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep2/Input_rep2_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed > Input_2_JF1_reads.sort.bed


~/tools_av/bedToBigBed Input_1_B6_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_1_B6_reads.bb
~/tools_av/bedToBigBed Input_2_B6_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_2_B6_reads.bb
~/tools_av/bedToBigBed Input_1_JF1_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_1_JF1_reads.bb
~/tools_av/bedToBigBed Input_2_JF1_reads.sort.bed /home/ankitv/ref_av/mm10/mm10.chrom.sizes Input_2_JF1_reads.bb

scp -r *.bb ankits@140.164.60.71:/media/ankits/Archivio_2/ankit/2019/chip-seq/Bigwigs/2019_new_chip-seq_peaks/
  
  
track type=bigBed name="JB1_Input_1_B6_reads" description="JB1_Input_1_B6_reads" bigDataUrl=http://140.164.60.71:8080/2019_new_chip-seq_peaks/Input_1_B6_reads.bb visibility=dense color=0,0,0
track type=bigBed name="JB1_Input_2_B6_reads" description="JB1_Input_2_B6_reads" bigDataUrl=http://140.164.60.71:8080/2019_new_chip-seq_peaks/Input_2_B6_reads.bb visibility=dense color=0,0,0
track type=bigBed name="JB1_Input_1_JF1_reads" description="JB1_Input_1_JF1_reads" bigDataUrl=http://140.164.60.71:8080/2019_new_chip-seq_peaks/Input_1_JF1_reads.bb visibility=dense color=0,0,0
track type=bigBed name="JB1_Input_2_JF1_reads" description="JB1_Input_2_JF1_reads" bigDataUrl=http://140.164.60.71:8080/2019_new_chip-seq_peaks/Input_2_JF1_reads.bb visibility=dense color=0,0,0

#Count Peaks
bedtools getfasta -fi /home/ankitv/ref_av/mm10_bowtie2/mm10.fa -bed NRF1_intersection_rep1_2_peaks.bed -fo NRF1_intersection_rep1_2_peaks.fa
/home/ankitv/meme/bin/meme-chip -meme-nmotifs 3 NRF1_intersection_rep1_2_peaks.fa

#ChIP-Peaks heatmap
#Obtain reference transcription start site TSS
http://reftss.clst.riken.jp/reftss/Main_Page
See description http://reftss.clst.riken.jp/datafiles/current/00README.txt

#Obtain reference transcription start site TSS
http://reftss.clst.riken.jp/reftss/Main_Page
See description http://reftss.clst.riken.jp/datafiles/current/00README.txt

#merge biwig of rep1 and rep2 #inside bigwigs folder
~/tools_av/bigWigMerge NRF1_treat_rep1_R1_mm10.bw NRF1_treat_rep2_R1_mm10.bw NRF1_treat_bulkmerge_R1_mm10.bdg -max
sort -k1,1 -k2,2n  NRF1_treat_bulkmerge_R1_mm10.bdg > NRF1_treat_bulkmerge_R1_mm10.sort.bdg
~/tools_av/bedGraphToBigWig NRF1_treat_bulkmerge_R1_mm10.sort.bdg /home/ankitv/ref_av/mm10/mm10.chrom.M.sizes  NRF1_treat_bulkmerge_R1_mm10.bw
computeMatrix reference-point -S NRF1_treat_bulkmerge_R1_mm10.bw -R refTSS_v3.1_mouse_coordinate.mm10.bed -o NRF1_treat_bulkmerge_R1_mm10_matrix.gz -a 3000 -b 3000 -bs 25 -p 15 --missingDataAsZero
plotHeatmap -m NRF1_treat_bulkmerge_R1_mm10_matrix.gz --colorList "white,green" -out NRF1_treat_bulkmerge_R1_mm10_matrix.png --sortUsing max

############################################   END OF ANALYSIS   ########################################################
bedtools coverage -a NRF1_intersection_rep1_2_peaks.bed -b JB1resultInputBi_chr.txt  > NRF1_intersection_rep1_2_peaks_covSNP.bed
bedtools intersect -wa -wb -a NRF1_intersection_rep1_2_peaks_covSNP.bed -b AlspNRF1ddsNormcounts_B6_monofilt_uniqB6.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""B6"}' > NRF1_intersection_rep1_2_peaks_covSNP_B6.bed
bedtools intersect -wa -wb -a NRF1_intersection_rep1_2_peaks_covSNP.bed -b AlspNRF1ddsNormcounts_JF1_monofilt_uniqJF1.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""JF1"}' > NRF1_intersection_rep1_2_peaks_covSNP_JF1.bed
cat NRF1_intersection_rep1_2_peaks_covSNP_B6.bed NRF1_intersection_rep1_2_peaks_covSNP_JF1.bed > NRF1_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed
bedtools intersect -wa -wb -a NRF1_intersection_rep1_2_peaks_covSNP.bed -b NRF1_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed -v | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t""Biallelic"}' > NRF1_intersection_rep1_2_peaks_covSNP_Biallelic.bed
cat NRF1_intersection_rep1_2_peaks_covSNP_B6_or_JF1.bed NRF1_intersection_rep1_2_peaks_covSNP_Biallelic.bed > NRF1_intersection_rep1_2_peaks_covSNP_B6_or_JF1_or_Biallelic.bed

JB1_NRF1_peaks_Het_SNP <- read.table("NRF1_intersection_rep1_2_peaks_covSNP_B6_or_JF1_or_Biallelic.bed", sep = "\t", header = F, stringsAsFactors = F)
JB1_NRF1_peaks_Het_SNP <- JB1_NRF1_peaks_Het_SNP[,c(1,2,3,6,7)]
head(JB1_NRF1_peaks_Het_SNP)
colnames(JB1_NRF1_peaks_Het_SNP) <- c("chr","start","end","SNP","allele")
head(JB1_NRF1_peaks_Het_SNP)
JB1_NRF1_peaks_Het_SNP["SNP_allele"] <- paste0(JB1_NRF1_peaks_Het_SNP$SNP, "_" ,JB1_NRF1_peaks_Het_SNP$allele)
countJB1_NRF1_peaks_Het_SNP <- count(JB1_NRF1_peaks_Het_SNP,"SNP_allele")
countJB1_NRF1_peaks_Het_SNPsplit <- cSplit(countJB1_NRF1_peaks_Het_SNP,"SNP_allele","_")
countJB1_NRF1_peaks_Het_SNPsplit <- countJB1_NRF1_peaks_Het_SNPsplit[order(countJB1_NRF1_peaks_Het_SNPsplit$SNP_allele_2),]

ggplot(data=countJB1_NRF1_peaks_Het_SNPsplit, aes(x=SNP_allele_1, y=freq)) +
  geom_bar(stat="identity", aes(fill=as.character(SNP_allele_2)),position=position_dodge())+
  theme_bw()+
  scale_fill_manual(values = c("darkred","darkgreen","navy"))
ggsave("countJB1_NRF1_peaks_Het_SNP.png", width=17*1.25, height=10*1.25, units="cm", dpi=96)

head(JB1result[["JB1"]])
SNPshetJB1result <- JB1result[["JB1"]]
rownames(SNPshetJB1result) <- paste0(SNPshetJB1result$CHROM, "_", SNPshetJB1result$POS)
SNPshetJB1result <- SNPshetJB1result[,c(6,7)]
SNPshetJB1resultst <- stack(as.matrix(SNPshetJB1result))
SNPshetJB1resultst <- data.frame(SNPshetJB1resultst)
ggplot(SNPshetJB1resultst, aes(x=value, fill=col, color=col)) +
  geom_histogram(position="dodge",alpha=0.5)+theme_bw()
ggsave("SNPshetJB1resultst.png", width=17*1.25, height=10*1.25, units="cm", dpi=96)

#Count SNPs
NRF1_B6_specific_peaks_jf1v2 <- read.table("NRF1_B6_specific_peaks_jf1v2.bed")
Count_NRF1_B6_jf1v2 <- count(NRF1_B6_specific_peaks_jf1v2, "V4")
head(Count_NRF1_B6_jf1v2)
dim(Count_NRF1_B6_jf1v2)
write.table(Count_NRF1_B6_jf1v2, "Count_NRF1_B6_jf1v2.txt", sep="\t", quote=F, append = F, row.names = F)

NRF1_JF1_specific_peaks_jf1v2 <- read.table("NRF1_JF1_specific_peaks_jf1v2.bed")
Count_NRF1_JF1_jf1v2 <- count(NRF1_JF1_specific_peaks_jf1v2, "V4")
head(Count_NRF1_JF1_jf1v2)
dim(Count_NRF1_JF1_jf1v2)
write.table(Count_NRF1_JF1_jf1v2, "Count_NRF1_JF1_jf1v2.txt", sep="\t", quote=F, append = F, row.names = F)


#Make a heatmap with all filtered peaks
AlspInputddsNormcounts_combinefilt <- rbind.data.frame(AlspInputddsNormcounts_JF1_monofilt, AlspInputddsNormcounts_B6_monofilt)

head(AlspInputddsNormcounts_combinefilt)

AlspInputddsNormcounts_combinefiltavg <- data.frame(cbind(data.frame(rownames(AlspInputddsNormcounts_combinefilt)),
                                                          (AlspInputddsNormcounts_combinefilt$chr),
                                                          (AlspInputddsNormcounts_combinefilt$start),
                                                          (AlspInputddsNormcounts_combinefilt$end),
                                                          (AlspInputddsNormcounts_combinefilt$Input_B6_1 + AlspInputddsNormcounts_combinefilt$Input_B6_2)/2,
                                                          (AlspInputddsNormcounts_combinefilt$Input_JF1_1 + AlspInputddsNormcounts_combinefilt$Input_JF1_2)/2))

head(AlspInputddsNormcounts_combinefiltavg)

rownames(AlspInputddsNormcounts_combinefiltavg) <- AlspInputddsNormcounts_combinefilt$peakID

AlspInputddsNormcounts_combinefiltavg <- AlspInputddsNormcounts_combinefiltavg[,-1]
head(AlspInputddsNormcounts_combinefiltavg)
colnames(AlspInputddsNormcounts_combinefiltavg) <- c("chr","start","end","InputB6","InputJF1")
head(AlspInputddsNormcounts_combinefiltavg)
AlspInputddsNormcounts_combinefiltavg["matRatioInput"] <- AlspInputddsNormcounts_combinefiltavg$InputJF1/(AlspInputddsNormcounts_combinefiltavg$InputB6 + AlspInputddsNormcounts_combinefiltavg$InputJF1)
head(AlspInputddsNormcounts_combinefiltavg)
dim(AlspInputddsNormcounts_combinefiltavg)

AlspInputddsNormcounts_combinefiltavgbins <- AlspInputddsNormcounts_combinefiltavg
head(AlspInputddsNormcounts_combinefiltavgbins)
dim(AlspInputddsNormcounts_combinefiltavgbins)
write.table(AlspInputddsNormcounts_combinefiltavgbins, "AlspInputddsNormcounts_combinefiltavgbins.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

library(circlize)
library(gtrellis)
AlspInputddsNormcounts_combinefiltavgbinsbed = AlspInputddsNormcounts_combinefiltavgbins[,c(1,2,3,6)]
gtrellis_layout(track_ylim = range(AlspInputddsNormcounts_combinefiltavgbinsbed[[4]]), nrow = 3, byrow = FALSE)
add_points_track(AlspInputddsNormcounts_combinefiltavgbinsbed, AlspInputddsNormcounts_combinefiltavgbinsbed[[4]], gp = gpar(col = ifelse(AlspInputddsNormcounts_combinefiltavgbinsbed[[4]] > 0, "red", "green")))

library(ComplexHeatmap)

lgd = Legend(at = c("class1", "class2"), title = "Class", type = "points", legend_gp = gpar(col = 2:3))
gtrellis_layout(nrow = 5, byrow = FALSE, track_ylim = range(bed[[4]]), legend = lgd)
add_points_track(bed, bed[[4]], gp = gpar(col = sample(2:3, nrow(bed), replace = TRUE)))
svg(filename="TABInputNRF1peak1.NRF1peak.svg", width=15, height=40, pointsize=12)

#Statistical test for allele specific
