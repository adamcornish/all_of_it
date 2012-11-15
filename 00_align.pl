#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopt ( "c:", \%opt );
my $config_file = $opt{c};
my $config      = `cat $config_file`;
my $ref_dir     = $config =~ /REF_DIR\s+(\S+)/;
my $index_dir   = $config =~ /INDEX_DIR\s+(\S+)/;
my $sample_type = $config =~ /SAMPLE_TYPE\s+(\S+)/; # can be "exome" or "rna-seq"
my $cancer      = $config =~ /CANCER\s+(\S+)/;      # can be "yes" or "no"
my $aligner     = $config =~ /ALIGNER\s+(\S+)/;     # can be "bwa", "bowtie2", or "both"
my $fasta       = $config =~ /FASTA\s+(\S+)/;
my $threads     = $config =~ /THREADS\s+(\S+)/;
my $reads_dir   = $config =~ /READS_DIR\s+(\S+)/;
my $exp_name    = $config =~ /NAME\s+(\S+)/;

chomp ( my @reads = `ls $reads_dir/*fastq` );

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name) = $reads[$i] =~ /^.+\/(.+?)(?:_|\.)/;
    my ($R1)   = $reads[$i] =~ /^.+\/(.+)/;
    my ($R2)   = $reads[$i+1] =~ /^.+\/(.+)/;
    if ( $aligner =~ /(?:bowtie2|both)/i )
    {
        chomp ( my ($btx) = `ls $index_dir/*.4.bt2 | sed s/\.4\.bt2//` );
        my $aln_method = ( $sample_type =~ /"exome"/i ) ? "--very-sensitive" : "--very-sensitive-local";
        my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bt2
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=8GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.stdout
cd $reads_dir
bowtie2 --rg-id $name --rg PL:illumina --rg PU:$name --rg LB:$name --rg SM:$name $aln_method -x $btx -p $threads -1 $R1 -2 $R2 -S tmp/$name.bt2.sam
perl 01_sam_to_bam.pl -c $config_file
END
        open OUT, ">qsub/00_$name.bt2_aln.qsub";
        print OUT $sub_script;
        close OUT;
    }
    if ( $aligner =~ /(?:bwa|both)/i )
    {
        my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=8GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.stdout
cd $reads_dir
bwa aln $fasta $R1 > tmp/$name.r1.sai
bwa aln $fasta $R2 > tmp/$name.r2.sai
bwa sampe -r '\@RG\\tID:$name\\tPL:illumina\\tPU:$name\\tLB:$name\\tSM:$name' $fasta tmp/$name.r1.sai tmp/$name.r2.sai $R1 $R2 > tmp/$name.bwa.sam
perl 01_sam_to_bam.pl -c $config_file
END
        open OUT, ">qsub/00_$name.bwa_aln.qsub";
        print OUT $sub_script;
        close OUT;
    }
   #system ( "qsub *$name*_aln.qsub*" ); # doing this will make it so we submit either one or two qsub scripts
}

#TODO add tophat for rnaseq to do fusion detection
