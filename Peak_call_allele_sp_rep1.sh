#Peak calling: (Allele separate)
#Concatenate Input B6 + JF1
cat Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed Input_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed >  Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6_JF1.rmdup_chr.bed 

#Peak calling: B6+JF1 (Allele merged)
cd NFY
cd rep1
cd merge_input_PeakCall_MACS1
macs14 -t ./../NFY_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6_JF1.rmdup_chr.bed -g mm -n NFY_rep1_B6 -B -S --call-subpeaks --keep-dup=all 

macs14 -t ./../NFY_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6_JF1.rmdup_chr.bed -g mm -n NFY_rep1_JF1 -B -S --call-subpeaks --keep-dup=all 

cd ..
cd ..
cd ..

cd NRF1
cd rep1
cd merge_input_PeakCall_MACS1/
macs14 -t ./../NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.B6.rmdup_chr.bed -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6_JF1.rmdup_chr.bed -g mm -n NRF1_rep1_B6 -B -S --call-subpeaks --keep-dup=all 

macs14 -t ./../NRF1_rep1_R1_uniq.sortedByReadname.genome.sort.JF1.rmdup_chr.bed -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/allele_sp/Input/rep1/Input_rep1_R1_uniq.sortedByReadname.genome.sort.B6_JF1.rmdup_chr.bed -g mm -n NRF1_rep1_JF1 -B -S --call-subpeaks --keep-dup=all 

cd ..
cd ..
cd ..

echo "Peak calling finished for rep1 Allele specific data"


