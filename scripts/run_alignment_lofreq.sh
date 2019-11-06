#!/bin/bash

#List of used software:
# - Trimmomatic
# - BWA
# - Picard Tools
# - GATK
# - Samtools
# - lofreq
# - HTSlib (for bgzip and tabix)
# - VCFtools
# - SnpEff

#The variant calling steps:
# -1- Install all software from the list
# -2- Create Ecoli_BW25113 database for SnpEff using included Escherichia_coli_BW25113.gbff GenBank file
# -3- Run prepare_reference.sh to create all index files for the reference
# -4- Copy fastq.gz files from the SRA
# -5- Create working directory
# -6- Copy run_alignment_lofreq.sh to the working directory and run it
# -7- Run the run_alignment_lofreq.sh

#specify according to your system
trimmdir=""			#path to Trimmomatic
bwadir=""			#path to BWA
picarddir=""		#path to Picard Tools
gatkdir=""			#path to GATK
vcfdir=""			#path to vcftools
snpeffdir=""		#path to snpEff
aux=""				#path to the auxiliary file folder (Aux in the archive)
fqdir=""			#directory with FastQ files
thr=1				#number of CPU threads

#create soft links for fastq files
for fq in $(ls $fqdir)
do
	dir=${fq%%_*}
	if [ ! -d $dir ]
	then
		mkdir $dir
	fi
	
	ln -s $fqdir/$fq $dir/$fq
done

#run trimming process

for d in $(ls -d */ | grep '[0-9]')
do
cd ./$d
    for r1 in $(ls *R1_001.fastq.gz)
    do
    	r1p=${r1/.fastq.gz/_1.fastq.gz}				#trimmed reads, 1st in pair
    	r2=${r1/R1/R2}
    	r2p=${r2/.fastq.gz/_1.fastq.gz}				#trimmed reads, 2nd in pair
    	r1s=${r1/.fastq.gz/_single_1.fastq.gz}		#trimmed reads, 1st in pair, pair absent
    	r2s=${r2/.fastq.gz/_single_1.fastq.gz}		#trimmed reads, 1st in pair, pair absent

    	java -jar $trimmdir/trimmomatic-0.36.jar PE -threads $thr -phred33 $r1 $r2 $r1p $r1s $r2p $r2s ILLUMINACLIP:$aux/overrepresented.fasta:2:30:10 SLIDINGWINDOW:4:15 HEADCROP:10 CROP:145 MINLEN:65
    done
cd ..
done

#run alignment and variant calling

for d in $(ls -d */ | grep '[0-9]')
do
    cd ./$d
    cp $aux/make_RG_bwa.pl .		#copy script that generates read group information BWA
    cp $aux/known.vcf .				#copy known variants (observed in unevolved) for the base quality recalibration
	
	#BWA alignment
	for r1 in *R1_001_1.fastq.gz
	do
	    j=${r1/.fastq.gz/}
	    rg=$(perl make_RG_bwa.pl "$r1")
	    bwa_sam="$j.sam"
	    bwa_bam="$j.bam"
	    r2=${r1/R1/R2}
	
	    $bwadir/bwa mem -t $thr -M -R "$rg" $aux/Escherichia_coli_BW25113.fna $r1 $r2 > $bwa_sam
	    java -jar $picarddir/picard.jar SortSam I=$bwa_sam O=$bwa_bam SO=coordinate CREATE_INDEX=true
	    rm $bwa_sam
	done
	
	#Merge aligned BAM files
	samtools merge Sample.bam *R1_001_1.bam
	java -jar $picarddir/picard.jar SortSam I=Sample.bam O=Sample_i.bam SO=coordinate CREATE_INDEX=true
	
	if [ -e Sample_i.bam ]
	then
	    rm Sample.ba*
	    rm AO-*.ba*
	fi
	
	#realign reads near indels
	lofreq viterbi --ref $aux/Escherichia_coli_BW25113.fna --keepflags -o Sample_v.bam Sample_i.bam
	samtools sort -o Sample_vi.bam Sample_v.bam
	samtools index -b Sample_vi.bam
	
	if [ -e Sample_vi.bam ]
	then
	    rm Sample_i.ba*
	    rm Sample_v.ba*
	    rm Sample.ba*
	fi
	
	#reacalibrate base scores
	java -jar $gatkdir/GenomeAnalysisTK.jar -T BaseRecalibrator -R $aux/Escherichia_coli_BW25113.fna -I Sample_vi.bam -knownSites known.vcf -o Sample_re.bqsr.grp
	java -jar $gatkdir/GenomeAnalysisTK.jar -T PrintReads -R $aux/Escherichia_coli_BW25113.fna -I Sample_vi.bam -BQSR Sample_re.bqsr.grp -o Sample_r_bqsr.bam
	
	if [ -e Sample_r_bqsr.bam ]
	then
		rm Sample_re.bqsr.grp
	    rm Sample_vi.ba*
	    rm Sample_re.ba*
	    rm target_intervals.list
	fi
	
	#lofreq variant calling
	lofreq call-parallel --pp-threads $thr --call-indels -f $aux/Escherichia_coli_BW25113.fna -o variants.lofreq.vcf Sample_r_bqsr.bam
    cd ..
done

if [ ! -d "Results" ];
then
    mkdir Results
fi

#copy results
for d in $(ls -d */ | grep '[0-9]')
do
    cd ./$d
    i=${d%'/'}
    vcf="$d.vcf"
    cp ./variants.lofreq.vcf ../Results/$vcf
    cd ..
done

#annotate results
cd ./Results
cp $aux/ann_pars.pl .		#copy merged annotation vcf file prser
cp $aux/Escherichia_coli_BW25113.gff .		#copy GFF file with the genome annotation
cp $aux/r_coordinates.txt .	#copy file with the list of repeated region
cp $aux/gff.py .	#copy a GFF library
cp $aux/timeline.py .	#copy script that builds table

filelist=""

for file in *.vcf
do
    gz=${file/.vcf/.gz}
    echo "$file $gz"
    bgzip -c $file > $gz
    tabix -p vcf $gz
    filelist="$filelist $gz"
done

export PERL5LIB=$PERL5LIB:$vcfdir/vcftools_0.1.13/perl
vcf-merge $filelist > merged.gz
bgzip -d merged.gz
rm *.gz
java -jar $snpeffdir/snpEff/snpEff.jar -ud 100 -v Ecoli_BW25113 merged > ann_1.vcf
rm *.tbi
perl ann_pars.pl 
python3 timeline.py --gff Escherichia_coli_BW25113.gff --ann ann.vcf --rep r_coordinates.txt --vcf_folder --out results.txt