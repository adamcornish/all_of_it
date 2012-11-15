#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

## Sort bam
# sorts the 

my %opt;
getopt ( "c:n:", \%opt );
my $name        = $opt{n};
my $config      = `cat $opt{c}`;
my ($bin)         = $config =~ /BIN\s+(\S+)/;
my ($index_dir)   = $config =~ /INDEX_DIR\s+(\S+)/;
my ($annotator)   = $config =~ /ANNOTATOR\s+(\S+)/;   # can be "snpeff" or "annovar"
my ($sample_type) = $config =~ /SAMPLE_TYPE\s+(\S+)/; # can be "exome" or "rna-seq"
my ($cancer)      = $config =~ /CANCER\s+(\S+)/;      # can be "yes", "no", "true", or "false"
my ($aligner)     = $config =~ /ALIGNER\s+(\S+)/;     # can be "bwa", "bowtie2", or "both"
my ($genotyper)   = $config =~ /GENOTYPER\s+(\S+)/;   # can be "gatk", "mpileup", or "both"
my ($fasta)       = $config =~ /FASTA\s+(\S+)/;
my ($threads)     = $config =~ /THREADS\s+(\S+)/;
my ($reads_dir)   = $config =~ /READS_DIR\s+(\S+)/;
my ($exp_name)    = $config =~ /NAME\s+(\S+)/;
my ($dbsnp)       = $config =~ /DBSNP\s+(\S+)/;
my ($omni)        = $config =~ /OMNI\s+(\S+)/;
my ($hapmap)      = $config =~ /HAPMAP\s+(\S+)/;
my ($mills)       = $config =~ /MILLS\s+(\S+)/;
my $gatk          = "$bin/GenomeAnalysisTK.jar";

if ( $aligner =~ /(?:bowtie2|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bt2.sam_to_bam
##PBS -l select=1
##PBS -l mem=4GB
##PBS -l walltime=2:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.sam_to_bam.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.sam_to_bam.stdout
cd $reads_dir
samtools view -bS tmp/$name.bt2.sam -o tmp/$name.bt2.bam
perl 02_sort_sam.pl -c $config_file -n $name
END
    open OUT, ">qsub/01_$name.bt2.sam_to_bam.qsub";
    print OUT $sub_script;
    close OUT;
}
if ( $aligner =~ /(?:bwa|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa.sam_to_bam
##PBS -l select=1
##PBS -l mem=4GB
##PBS -l walltime=2:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.sam_to_bam.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.sam_to_bam.stdout
cd $reads_dir
samtools view -bS tmp/$name.bwa.sam -o tmp/$name.bwa.bam
perl 02_sort_sam.pl -c $config_file -n $name
END
    open OUT, ">qsub/01_$name.bwa.sam_to_bam.qsub";
    print OUT $sub_script;
    close OUT;
}
#system ( "qsub *$name*_aln.qsub*" ); # doing this will make it so we submit either one or two qsub scripts
