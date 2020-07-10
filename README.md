# Triclosan Experimental evolution paper
Experimental evolution of E. coli with exposure to triclosan

The repository contains code that was used in the data analysis for experimental evolution of *Escherichia coli BW25113 &Delta;uxaC::kan* strain under the triclosan treatment.

### List of used software:
- Trimmomatic ([site](http://www.usadellab.org/cms/?page=trimmomatic))
- BWA ([site](http://bio-bwa.sourceforge.net/))
- Picard Tools ([site](https://broadinstitute.github.io/picard/))
- GATK ([site](https://software.broadinstitute.org/gatk/))
- Samtools and HTSlib ([site](http://www.htslib.org/))
- lofreq ([site](https://csb5.github.io/lofreq/))
- VCFtools ([site](http://vcftools.sourceforge.net/))
- SnpEff ([site](http://snpeff.sourceforge.net/))

### The variant calling steps:
1. Install all software from the list.
2. Create Ecoli_BW25113 database for SnpEff using included `Escherichia_coli_BW25113.gbff` GenBank file (should be renamed to `genes.gbk`).
3. Run `prepare_reference.sh` to create all index files for the reference.
4. Copy fastq.gz files from the SRA.
5. Create working directory.
6. Copy `run_alignment_lofreq.sh` to the working directory.
7. Fill all required paths in `run_alignment_lofreq.sh`.
8. Run the `run_alignment_lofreq.sh`.
