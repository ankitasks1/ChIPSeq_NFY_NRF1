echo "Peak calling: NFY rep2"
cd rep2
mkdir PeakCall_MACS1
cd PeakCall_MACS1
macs14 -t ./../NFY_rep2_R1_uniq.rmdup.bam  -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/Input/Input_rep2_R1_uniq.rmdup.bam -g mm -n NFY_rep2 -B -S --call-subpeaks --keep-dup all
cd ..
mkdir PeakCall_MACS2
cd PeakCall_MACS2
macs2 callpeak -t ./../NFY_rep2_R1_uniq.rmdup.bam -c /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bulk/Input/Input_rep2_R1_uniq.rmdup.bam -g mm -n NFY_rep2 -B --call-summits --keep-dup all
cd ..
cd ..







