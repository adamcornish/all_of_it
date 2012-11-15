#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopt ( "c:", \%opt );
my $config       = `cat $opt{c}`;
my $bin          = ( $config =~ /BIN\s+(\S+)/ )         ? $1 : "/opt/var_calling";
my $ref_dir      = ( $config =~ /REF_DIR\s+(\S+)/ )     ? $1 : "/safer/genomes/Homo_sapiens/UCSC/hg19";
my $index_dir    = ( $config =~ /INDEX_DIR\s+(\S+)/ )   ? $1 : "$ref_dir/Sequence/BowtieIndex";
my $sample_type  = ( $config =~ /SAMPLE_TYPE\s+(\S+)/ ) ? $1 : "exome";   # can be "exome" or "rna-seq"
my $cancer       = ( $config =~ /CANCER\s+(\S+)/ )      ? $1 : "no";      # can be "yes" or "no"
my $aligner      = ( $config =~ /ALIGNER\s+(\S+)/ )     ? $1 : "bowtie2"; # can be "bwa", "bowtie2", or "both"
my $fasta        = ( $config =~ /FASTA\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
my $threads      = ( $config =~ /THREADS\s+(\S+)/ )     ? $1 : "24";
my $reads_dir    = ( $config =~ /READS_DIR\s+(\S+)/ )   ? $1 : ".";
my $exp_name     = ( $config =~ /NAME\s+(\S+)/ )        ? $1 : "serenity";

if ( $aligner eq "bowtie" )
{

    system ( "bowtie2 --very-sensitive-local -x $index_dir/-p $threads -1 $R1 -2 $R2 -S tmp/$name.sam" );
}
elsif ( $aligner eq "bwa" )
{
}
else
{
}
