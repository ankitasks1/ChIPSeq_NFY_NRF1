echo "Allele-Specific Analysis"
echo "Analysis pipeline: Single-end Mode Input"
#~/tools_av/FastQC/fastqc *.fastq.gz

echo "Trimming the original data"
#java -jar /home/ankitv/tools_av/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 Input_rep1_R1.fastq.gz Input_rep1_R1_trimmed.fastq.gz ILLUMINACLIP:/home/ankitv/tools_av/Trimmomatic-0.36/adapters/Combine-SE.fa:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#~/tools_av/FastQC/fastqc *.fastq.gz

echo "Aligning Data"
(/home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -p 10 -x /home/ankitv/ref_av/Nmm10_b2/Nmm10 -U Input_rep1_R1_trimmed.fastq.gz -S Input_rep1_R1.sam)2>Input_rep1_R1.log

echo "SAM to BAM"
samtools view -S -b Input_rep1_R1.sam > Input_rep1_R1.bam

echo "Extracting Uniquely Mapped Reads"
samtools view -b -q 20 Input_rep1_R1.bam > Input_rep1_R1_uniq.bam

echo "Sort BAM file by readname"
samtools sort -n Input_rep1_R1_uniq.bam -o Input_rep1_R1_uniq.sortedByReadname.bam 

echo "Splitting reads"
~/tools_av/SNPsplit_v0.3.2/SNPsplit --snp_file /home/ankitv/ref_av/jf1_ref/SNPsplit/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt Input_rep1_R1_uniq.sortedByReadname.bam

echo "Sort strain BAM files"
samtools sort -o Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.bam Input_rep1_R1_uniq.sortedByReadname.genome1.bam
samtools sort -o Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.bam Input_rep1_R1_uniq.sortedByReadname.genome2.bam


echo "Removing duplicates using picard"
java -jar /home/ankitv/tools_av/picard-2.26.2/picard.jar MarkDuplicates I=Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.bam O=Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam M=Input_rep1_R1_B6_dup.txt REMOVE_DUPLICATES=true CREATE_INDEX=true
java -jar /home/ankitv/tools_av/picard-2.26.2/picard.jar MarkDuplicates I=Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.bam O=Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup.bam M=Input_rep1_R1_JF1_dup.txt REMOVE_DUPLICATES=true CREATE_INDEX=true

echo "Convert BAM to Bed"
bedtools bamtobed -i Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.bam > Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.bed
bedtools bamtobed -i Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.bam > Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.bed

echo "Removing duplicates using bed"
cp Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.bed 1B6.bed
awk '{print $1"%"$2"%"$3"%"$6"\t"$4"\t"$5}' 1B6.bed > 2B6.bed
rm -r 1B6.bed
sort -k1,1 -k3,3rn 2B6.bed > 3B6.bed
rm -r 2B6.bed
sort -k1,1 -u 3B6.bed > 4B6.bed
rm -r 3B6.bed
awk -F "%" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' 4B6.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' > 5B6.bed
rm -r 4B6.bed
mv 5B6.bed Input_rep1_R1_uniq.sortedByReadname.genome.sort_B6_rmdup.bed

cp Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.bed 1JF1.bed
awk '{print $1"%"$2"%"$3"%"$6"\t"$4"\t"$5}' 1JF1.bed > 2JF1.bed
rm -r 1JF1.bed
sort -k1,1 -k3,3rn 2JF1.bed > 3JF1.bed
rm -r 2JF1.bed
sort -k1,1 -u 3JF1.bed > 4JF1.bed
rm -r 3JF1.bed
awk -F "%" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' 4JF1.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' > 5JF1.bed
rm -r 4JF1.bed
mv 5JF1.bed Input_rep1_R1_uniq.sortedByReadname.genome.sort_JF1_rmdup.bed

rm *.sam

echo "Analysis finished for Input"
