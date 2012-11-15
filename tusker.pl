#!/usr/bin/perl
use warnings;
use strict;

usage() unless $#ARGV == 0;
my $in     = shift;
my $config = `cat $in`;
usage() if $config !~ /READS_DIR/;

######## Start Variables ########

my $bin          = ( $config =~ /BIN\s+(\S+)/ )         ? $1 : "/opt/var_calling";
my $ref_dir      = ( $config =~ /REF_DIR\s+(\S+)/ )     ? $1 : "/safer/genomes/Homo_sapiens/UCSC/hg19";
my $index_dir    = ( $config =~ /INDEX_DIR\s+(\S+)/ )   ? $1 : "$ref_dir/Sequence/BowtieIndex";
my $annotator    = ( $config =~ /ANNOTATOR\s+(\S+)/ )   ? $1 : "snpeff";  # can be "snpeff" or "annovar"
my $sample_type  = ( $config =~ /SAMPLE_TYPE\s+(\S+)/ ) ? $1 : "exome";   # can be "exome" or "rna-seq"
my $cancer       = ( $config =~ /CANCER\s+(\S+)/ )      ? $1 : "no";      # can be "yes", "no", "true", or "false"
my $aligner      = ( $config =~ /ALIGNER\s+(\S+)/ )     ? $1 : "bowtie2"; # can be "bwa", "bowtie2", or "both"
my $genotyper    = ( $config =~ /GENOTYPER\s+(\S+)/ )   ? $1 : "gatk";    # can be "gatk", "mpileup", or "both"
my $fasta        = ( $config =~ /FASTA\s+(\S+)/ )       ? $1 : "$ref_dir/Sequence/WholeGenomeFasta/ucsc.hg19.fasta";
my $threads      = ( $config =~ /THREADS\s+(\S+)/ )     ? $1 : "24";
my $reads_dir    = ( $config =~ /READS_DIR\s+(\S+)/ )   ? $1 : ".";
my $exp_name     = ( $config =~ /NAME\s+(\S+)/ )        ? $1 : "serenity";
my $cluster      = ( $config =~ /CLUSTER\s+(\S+)/ )     ? $1 : "yes";     # can be "yes", "no", "true", or "false"
my $gatk         = "$bin/GenomeAnalysisTK.jar";
my $dbsnp        = "$ref_dir/Annotation/Variation/dbsnp.vcf";
my $omni         = "$ref_dir/Annotation/Variation/omni.vcf";
my $hapmap       = "$ref_dir/Annotation/Variation/hapmap.vcf";
my $mills        = "$ref_dir/Annotation/Variation/indels.vcf";
chomp ( my $time = `date +%T` );

run_sanity_checks ( $bin, $ref_dir, $index_dir, $annotator, $sample_type, $cancer, $aligner, $genotyper, $fasta, $threads, $reads_dir, $cluster );

print "Options used         :\n",
      "\tBIN          : $bin\n",
      "\tANNOTATOR    : $annotator\n",
      "\tSAMPLE_TYPE  : $sample_type\n",
      "\tCANCER       : $cancer\n",
      "\tALIGNER      : $aligner\n",
      "\tGENOTYPER    : $genotyper\n",
      "\tFASTA        : $fasta\n",
      "\tTHREADS      : $threads\n",
      "\tREADS_DIR    : $reads_dir\n",
      "\tNAME         : $exp_name\n";

system ( "00_run_aligner.pl -c $in" );

sub usage
{
    die <<USAGE;

    Usage: perl $0 <configuration_file.txt>

    Your configuration file MUST have a READS_DIR specified.

    Configuration options available

      OPTION      Default                                             Description
      BIN         /opt/var_calling                                    Absolute location of the Picard Tools and GATK jar files
      REF_DIR     /safer/genomes/Homo_sapiens/UCSC/hg19               Absolute location of the reference directory
      INDEX_DIR   REF_DIR/Sequence/BowtieIndex                        Absolute location of your reference indexes
      FASTA       REF_DIR/Sequence/WholeGenomeFasta/ucsc.hg19.fasta   Absolute location of the reference fasta file
      READS_DIR   .                                                   Absolute location of the reads that are going to be used
      ANNOTATOR   snpeff                                              Annotator to use: snpeff or annovar
      SAMPLE_TYPE exome                                               Sample type: exome or rna-seq
      CANCER      yes                                                 Whether or not this is a cancer sample: yes, no, true, false
      ALIGNER     bowtie2                                             Aligner to use: bowtie2, bwa, or both
      GENOTYPER   gatk                                                Genotyper to use: GATK, mpileup, or both
      THREADS     24                                                  Number of threads to use in parallelizable modules
      NAME        serenity                                            The name you want to give to this experiment

USAGE
}

sub run_sanity_checks
{
    my $bin         = shift;
    my $ref_dir     = shift;
    my $index_dir   = shift;
    my $annotator   = shift;
    my $sample_type = shift;
    my $cancer      = shift;
    my $aligner     = shift;
    my $genotyper   = shift;
    my $fasta       = shift;
    my $threads     = shift;
    my $reads_dir   = shift;
    my $cluster     = shift;
    chomp ( my $procs = `cat /proc/cpuinfo | grep -c processor` );

    die "The 'BIN' directory, $bin,  does not exist. Exiting.\n"                              unless ( -d $bin );
    die "The 'REF_DIR' directory, $ref_dir,  does not exist. Exiting.\n"                      unless ( -d $ref_dir );
    die "The 'INDEX_DIR' directory, $index_dir,  does not exist. Exiting.\n"                  unless ( -d $index_dir );
    die "Unknown annotator: your annotator options are 'annovar' or 'snpeff'. Exiting.\n"     unless ( $annotator =~ /(?:annovar|snpeff)/i );
    die "Unknown sample type: your sample type options are 'exome' or 'rna-seq'. Exiting.\n"  unless ( $sample_type =~ /(?:exome|rna-seq|rnaseq)/i );
    die "Unknown cancer option: your options are 'yes', 'no', 'true', or 'false'. Exiting.\n" unless ( $cancer =~ /(?:yes|no|true|false)/i );
    die "Unknown cluster option: your options are 'yes', 'no', 'true', or 'false'. Exiting.\n" unless ( $cancer =~ /(?:yes|no|true|false)/i );
    die "Unknown aligner: your optios are bwa, bowtie2, or both. Exiting.\n"                  unless ( $aligner =~ /(?:bwa|bowtie2|both)/i );
    die "Unknown genotyper: your options are GATK, mpileup, or both. Exiting.\n"              unless ( $genotyper =~ /(?:gatk|mpileup|both)/i );
    die "The genome fasta file, $fasta, does not exist. Exiting.\n"                           unless ( -e $fasta );
    die "The number of threads must be between 1 and $procs. Exiting.\n"                      unless ( $threads > 0 and $threads <= $procs );
    die "The 'READS_DIR' directory, $reads_dir, does not exist. Exiting.\n"                   unless ( -d $reads_dir );
    die "The 'READS_DIR' directory, $reads_dir, does not contain any fastq files. Exiting.\n" unless ( `ls $reads_dir/*fastq 2> /dev/null` );

    my ($ext) = $fasta =~ /(\.fasta|\.fa)/;
    my ($fa_name) = $fasta =~ /.+\/(.+?)$ext/;
    unless ( `ls $index_dir/$fa_name*bt2 2> /dev/null` )
    {
        system ( "mkdir logs" ) unless -d "logs";
        print "The 'INDEX_DIR' directory, $index_dir, does not contain any bowtie2 index files. Creating them now.\n";
        system ( "bowtie2-build $fasta $index_dir/$fa_name > logs/bowtie2-build.log 2> logs/bowtie2-build.log" );
    }
    unless ( `ls $index_dir/*bwt 2> /dev/null` )
    {
        system ( "mkdir logs" ) unless -d "logs";
        print "The 'INDEX_DIR' directory, $index_dir, does not contain any bwa index files. Creating them now.\n";
        system ( "bwa index -a bwtsw $fasta > logs/bwa-index.log 2> logs/bwa-index.log" );
    }
}
