#Effective genome size was used default: https://www.biostars.org/p/55271/
#==> samtools flagstat NFY_rep1_R1_uniq.rmdup.bam <== 13804282

#==> samtools flagstat NFY_rep2_R1_uniq.rmdup.bam <== 12165404

#==> samtools flagstat Input_rep1_R1_uniq.rmdup.bam <== 11139369

#==> samtools flagstat Input_rep2_R1_uniq.rmdup.bam <== 19590958

import os,sys

x = os.getcwd()

print x
print 'Analysis Started'
print '----------------'


lists ='''
rep1 Input Input NFY NFY 200,53,50 NFY 13804282 treat
rep1 Input Input NFY NFY 0,0,0 NFY 11139369 control
rep2 Input Input NFY NFY 200,53,50 NFY 12165404 treat
rep2 Input Input NFY NFY 0,0,0 NFY 19590958 control
'''.strip().split('\n')

for listx in lists:
    listx = listx.split(' ')
    folder = listx[0]
    gse_i = listx[1]
    id_i = listx[2]
    gse_s = listx[3]
    id_s = listx[4]
    color = listx[5]
    protein = listx[6]
    readcount = listx[7]
    peakfolder = listx[8]
    print folder
    print protein
    print color
    print 'Sample in analysis: ' + gse_s
    print 'Sample_SRR: '+ id_s
    print 'Sample_type: '+ peakfolder
    os.chdir(folder)
    #print os.listdir('./')
    header_bedgraph = 'track type=bedGraph name="'+ protein + '_' + gse_s + '_' + peakfolder + '_bdg_mm10" description="'+ protein + '_' + gse_s + '_' + peakfolder + '_bdg_mm10" visibility=full color=' + color
    #header_wig = 'track type=wiggle_0 name="' + protein + '_' + gse_s +'_wig_mm10" description="'+ protein + '_' + gse_s + '_wig_mm10" visibility=full color=' + color + '\n'
    #print header_bedgraph
    for sample in os.listdir('./'):
        if sample.endswith('_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam'):
           print 'Total number of reads in sample are: '
           os.system('samtools flagstat ' + protein + '_' + folder + '_R1_uniq.sortedByReadname.genome.sort.B6.rmdup.bam')
           print 'But original factor from Bulk is: ' + readcount + ', So use this'
           SF = (float(1000000)/float(readcount))
           print SF
           os.chdir('merge_input_PeakCall_MACS1')
           print 'Performing RPM normalisation'
           os.chdir(gse_s + '_' + folder + '_B6_MACS_bedGraph')
           os.chdir(peakfolder)
           os.system('cp *.bdg.gz peak.bdg.gz')
           os.system('gunzip peak.bdg.gz')
           os.system('sed -i -e "1d" peak.bdg')
           os.system('~/tools_av/macs_correction_mm10.sh')
           print 'Performing sample RPM normalisation'
           bedgraph_sample = open('corrected_peaks1.bdg','r')
           bedgraph_sample = bedgraph_sample.read().strip().split('\n')
           
           print 'Generating Bedgraph..'           
           with open(protein + '_' + peakfolder + '_' + folder + '_R1_B6_mm10.merge.bedgraph', 'a') as myfile:
                myfile.write(header_bedgraph + '\n')
           for lines in bedgraph_sample:
               lines = lines.split('\t')
               cov = str(int(lines[3]) * SF)
               my_bedgraph = (lines[0], lines[1], lines[2], cov)
               my_bedgraph = "\t".join(my_bedgraph)

               with open(protein + '_' + peakfolder + '_' + folder + '_R1_B6_mm10.merge.bedgraph', 'a') as myfile:
                    myfile.write(my_bedgraph +'\n')

           print 'Generating Bigwig...'
           os.system('cp ' + protein + '_' + peakfolder + '_' + folder + '_R1_B6_mm10.merge.bedgraph ' + protein + '_' + peakfolder + '_' + folder + '_B6_norm.bedgraph')
           os.system('sed -i -e "1d" *_B6_norm.bedgraph')
           os.system('sort -k1,1 -k2,2n ' + protein + '_' + peakfolder + '_' + folder + '_B6_norm.bedgraph > ' + protein + '_' + peakfolder + '_' + folder + '_B6_norm_sorted.bedgraph')
           os.system('~/tools_av/bedGraphToBigWig ' + protein + '_' + peakfolder + '_' + folder + '_B6_norm_sorted.bedgraph ' + '/home/ankitv/ref_av/mm10/mm10.chrom.sizes ' + protein + '_' + peakfolder + '_' + folder + '_R1_B6_mm10.merge.bw')
           os.system('gzip ' + protein + '_' + peakfolder + '_' + folder + '_R1_B6_mm10.merge.bedgraph')
           os.system('cp *.bw /media/ankitv/Archivio1/2019/chip-seq/JB1_NFY_NRF1/newway/bigwigs/')
           os.system('rm -r corrected_peaks1.bdg')
           os.system('rm -r *_B6_norm_sorted.bedgraph')
           os.system('rm -r *_B6_norm.bedgraph')
           print 'Analysis completed for ' + protein + '_' + folder + '_B6_' + peakfolder
           print '\n'         

    os.chdir('..')
    os.chdir('..')
    os.chdir('..')
    os.chdir('..')
