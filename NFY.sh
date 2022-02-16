cd rep1

echo "Analysis pipeline: Single-end Mode NFY rep1"
#~/tools_av/FastQC/fastqc *.fastq.gz

echo "Trimming the original data"
java -jar /home/ankitv/tools_av/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 NFY_rep1_R1.fastq.gz NFY_rep1_R1_trimmed.fastq.gz ILLUMINACLIP:/home/ankitv/tools_av/Trimmomatic-0.36/adapters/Combine-SE.fa:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#~/tools_av/FastQC/fastqc *_trimmed.fastq.gz

echo "Aligning Data"
(/home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -p 10 -x /home/ankitv/ref_av/mm10_bowtie2/mm10 -U NFY_rep1_R1_trimmed.fastq.gz -S NFY_rep1_R1.sam)2>NFY_rep1_R1.log
echo "SAM to BAM"
samtools view -S -b NFY_rep1_R1.sam > NFY_rep1_R1.bam
echo "Extracting Uniquely Mapped Reads"
samtools view -b -q 20 NFY_rep1_R1.bam > NFY_rep1_R1_uniq.bam
echo "Sort BAM file"
samtools sort -o NFY_rep1_R1_uniq.sort.bam NFY_rep1_R1_uniq.bam
echo "Removing duplicates using picard"
java -jar /home/ankitv/tools_av/picard-2.26.2/picard.jar MarkDuplicates I=NFY_rep1_R1_uniq.sort.bam O=NFY_rep1_R1_uniq.rmdup.bam M=NFY_rep1_R1_dup.txt REMOVE_DUPLICATES=true CREATE_INDEX=true
echo "Convert BAM to Bed"
bedtools bamtobed -i NFY_rep1_R1_uniq.sort.bam > NFY_rep1_R1_uniq.sort.bed
echo "Removing duplicates using bed"
cp NFY_rep1_R1_uniq.sort.bed 1ext1.bed
awk '{print $1"%"$2"%"$3"%"$6"\t"$4"\t"$5}' 1ext1.bed > 2ext1.bed
rm -r 1ext1.bed
sort -k1,1 -k3,3rn 2ext1.bed > 3ext1.bed
rm -r 2ext1.bed
sort -k1,1 -u 3ext1.bed > 4ext1.bed
rm -r 3ext1.bed
awk -F "%" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' 4ext1.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' > 5ext1.bed
rm -r 4ext1.bed
mv 5ext1.bed NFY_rep1_R1_uniq.sort_rmdup.bed
rm -r *.sam

cd ..
cd rep2
echo "Analysis pipeline: Single-end Mode NFY rep2"
#~/tools_av/FastQC/fastqc *.fastq.gz

echo "Trimming the original data"
java -jar /home/ankitv/tools_av/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 NFY_rep2_R1.fastq.gz NFY_rep2_R1_trimmed.fastq.gz ILLUMINACLIP:/home/ankitv/tools_av/Trimmomatic-0.36/adapters/Combine-SE.fa:2:30:10 AVGQUAL:20 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#~/tools_av/FastQC/fastqc *_trimmed.fastq.gz

echo "Aligning Data"
(/home/ankitv/tools_av/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -p 10 -x /home/ankitv/ref_av/mm10_bowtie2/mm10 -U NFY_rep2_R1_trimmed.fastq.gz -S NFY_rep2_R1.sam)2>NFY_rep2_R1.log
echo "SAM to BAM"
samtools view -S -b NFY_rep2_R1.sam > NFY_rep2_R1.bam
echo "Extracting Uniquely Mapped Reads"
samtools view -b -q 20 NFY_rep2_R1.bam > NFY_rep2_R1_uniq.bam
echo "Sort BAM file"
samtools sort -o NFY_rep2_R1_uniq.sort.bam NFY_rep2_R1_uniq.bam
echo "Removing duplicates using picard"
java -jar /home/ankitv/tools_av/picard-2.26.2/picard.jar MarkDuplicates I=NFY_rep2_R1_uniq.sort.bam O=NFY_rep2_R1_uniq.rmdup.bam M=NFY_rep2_R1_dup.txt REMOVE_DUPLICATES=true CREATE_INDEX=true
echo "Convert BAM to Bed"
bedtools bamtobed -i NFY_rep2_R1_uniq.sort.bam > NFY_rep2_R1_uniq.sort.bed
echo "Removing duplicates using bed"
cp NFY_rep2_R1_uniq.sort.bed 1ext1.bed
awk '{print $1"%"$2"%"$3"%"$6"\t"$4"\t"$5}' 1ext1.bed > 2ext1.bed
rm -r 1ext1.bed
sort -k1,1 -k3,3rn 2ext1.bed > 3ext1.bed
rm -r 2ext1.bed
sort -k1,1 -u 3ext1.bed > 4ext1.bed
rm -r 3ext1.bed
awk -F "%" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' 4ext1.bed | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$4}' > 5ext1.bed
rm -r 4ext1.bed
mv 5ext1.bed NFY_rep2_R1_uniq.sort_rmdup.bed
rm -r *.sam

cd ..


