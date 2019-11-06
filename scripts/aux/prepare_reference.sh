bwadir=""			#path to BWA
picarddir=""		#path to Picard Tools

fasta="Escherichia_coli_BW25113.fna"
dict=${fasta/.fna/.dict/}	
$bwadir/bwa index $fasta
java -jar $picarddir/picard.jar CreateSequenceDictionary R=$fasta O=$dict
samtools faidx $fasta