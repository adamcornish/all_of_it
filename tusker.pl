#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopt ( "hc:", \%opt );
usage() if ( $opt{h} or !$opt{c} );
my $config_file = $opt{c};
my $config      = `cat $config_file`;

######## Start Variables ########
# TODO add option for dry run; i.e. do not do system ( "qsub submit_script.qsub" ); unless $dry_run = false

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
my $script        = "";
chomp ( my $time  = `date +%T` );

######## Sanity Checks ########

chomp ( my $procs = `cat /proc/cpuinfo | grep -c processor` );

# Make sure the config file is formatted correctly
die "BIN variable was not found in the configuration file, $config_file. Exiting.\n"         unless $bin;
die "INDEX_DIR variable was not found in the configuration file, $config_file. Exiting.\n"   unless $index_dir;
die "ANNOTATOR variable was not found in the configuration file, $config_file. Exiting.\n"   unless $annotator;
die "SAMPLE_TYPE variable was not found in the configuration file, $config_file. Exiting.\n" unless $sample_type;
die "CANCER variable was not found in the configuration file, $config_file. Exiting.\n"      unless $cancer;
die "ALIGNER variable was not found in the configuration file, $config_file. Exiting.\n"     unless $aligner;
die "GENOTYPER variable was not found in the configuration file, $config_file. Exiting.\n"   unless $genotyper;
die "FASTA variable was not found in the configuration file, $config_file. Exiting.\n"       unless $fasta;
die "THREADS variable was not found in the configuration file, $config_file. Exiting.\n"     unless $threads;
die "READS_DIR variable was not found in the configuration file, $config_file. Exiting.\n"   unless $reads_dir;
die "EXP_NAME variable was not found in the configuration file, $config_file. Exiting.\n"    unless $exp_name;
die "DBSNP variable was not found in the configuration file, $config_file. Exiting.\n"       unless $dbsnp;
die "OMNI variable was not found in the configuration file, $config_file. Exiting.\n"        unless $omni;
die "HAPMAP variable was not found in the configuration file, $config_file. Exiting.\n"      unless $hapmap;
die "MILLS variable was not found in the configuration file, $config_file. Exiting.\n"       unless $mills;
# Check to make sure everything exists
die "The 'BIN' directory, $bin,  does not exist. Exiting.\n"                                 unless ( -d $bin );
die "The 'INDEX_DIR' directory, $index_dir,  does not exist. Exiting.\n"                     unless ( -d $index_dir );
die "The 'READS_DIR' directory, $reads_dir, does not exist. Exiting.\n"                      unless ( -d $reads_dir );
die "The genome fasta file, $fasta, does not exist. Exiting.\n"                              unless ( -e $fasta );
die "The dbsnp vcf file, $dbsnp, does not exist. Exiting.\n"                                 unless ( -e $dbsnp);
die "The omni vcf file, $omni, does not exist. Exiting.\n"                                   unless ( -e $omni);
die "The hapmap vcf file, $hapmap, does not exist. Exiting.\n"                               unless ( -e $hapmap);
die "The mills vcf file, $mills, does not exist. Exiting.\n"                                 unless ( -e $mills);
die "The 'READS_DIR' directory, $reads_dir, does not contain any fastq files. Exiting.\n"    unless ( `ls $reads_dir/*fastq 2> /dev/null` );
# Validate formatting
die "Unknown annotator: your annotator options are 'annovar' or 'snpeff'. Exiting.\n"        unless ( $annotator =~ /(?:annovar|snpeff)/i );
die "Unknown sample type: your sample type options are 'exome' or 'rna-seq'. Exiting.\n"     unless ( $sample_type =~ /(?:exome|rna-seq|rnaseq)/i );
die "Unknown cancer option: your options are 'yes', 'no', 'true', or 'false'. Exiting.\n"    unless ( $cancer =~ /(?:yes|no|true|false)/i );
die "Unknown aligner: your options are bwa, bowtie2, or both. Exiting.\n"                    unless ( $aligner =~ /(?:bwa|bowtie2|both)/i );
die "Unknown genotyper: your options are GATK, mpileup, or both. Exiting.\n"                 unless ( $genotyper =~ /(?:gatk|mpileup|both)/i );
die "The number of threads must be between 1 and $procs. Exiting.\n"                         unless ( $threads > 0 and $threads <= $procs );
# Make sure that the locations of things are does as absolute paths, not relative paths
die "The 'READS_DIR' directory, $reads_dir, should be presented as an absolute path, not a relative path. Exiting.\n"       unless ( $reads_dir =~ /^\// );
die "The 'BIN' directory, $bin, should be presented as an absolute path, not a relative path. Exiting.\n"                   unless ( $bin =~ /^\// );
die "The 'INDEX_DIR' directory, $index_dir, should be presented as an absolute path, not a relative path. Exiting.\n"       unless ( $index_dir =~ /^\// );
die "The location of the 'FASTA' file, $fasta, should be presented as an absolute path, not a relative path. Exiting.\n"    unless ( $fasta =~ /^\// );
die "The location of the 'DBSNP' file, $dbsnp, should be presented as an absolute path, not a relative path. Exiting.\n"    unless ( $dbsnp =~ /^\// );
die "The location of the 'OMNI' file, $omni, should be presented as an absolute path, not a relative path. Exiting.\n"      unless ( $omni =~ /^\// );
die "The location of the 'HAPMAP' file, $hapmap, should be presented as an absolute path, not a relative path. Exiting.\n"  unless ( $hapmap =~ /^\// );
die "The location of the 'MILLS' file, $mills, should be presented as an absolute path, not a relative path. Exiting.\n"    unless ( $mills =~ /^\// );

# Create needed directories
system ( "mkdir $reads_dir/tmp" )  unless -d "$reads_dir/tmp";
system ( "mkdir $reads_dir/logs" ) unless -d "$reads_dir/logs";
system ( "mkdir $reads_dir/qsub" ) unless -d "$reads_dir/qsub";
my ($ext)     = $fasta =~ /(\.fasta|\.fa)/;
my ($fa_name) = $fasta =~ /.+\/(.+?)$ext/;

# Make index files if they don't exist
# TODO make sure to change this so it's run on a job node and not the submit node
unless ( `ls $index_dir/$fa_name*bt2 2> /dev/null` )
{
    print "The 'INDEX_DIR' directory, $index_dir, does not contain any bowtie2 index files for the fasta file $fa_name.$ext. Creating them now.\n";
    system ( "bowtie2-build $fasta $index_dir/$fa_name > $reads_dir/logs/bowtie2-build.log 2> $reads_dir/logs/bowtie2-build.log" );
}
unless ( `ls $index_dir/$fa_name*bwt 2> /dev/null` )
{
    print "The 'INDEX_DIR' directory, $index_dir, does not contain any bwa index files for the fasta file $fa_name.$ext. Creating them now.\n";
    system ( "bwa index -a bwtsw $fasta > $reads_dir/logs/bwa-index.log 2> $reads_dir/logs/bwa-index.log" );
}

## All checks have passed, so start the pipeline

chomp ( my @reads = `ls $reads_dir/*fastq` );

for ( my $i = 0; $i < @reads; $i += 2 )
{
    my ($name) = $reads[$i] =~ /^.+\/(.+?)(?:_|\.)/;
    my ($R1)   = $reads[$i] =~ /^.+\/(.+)/;
    my ($R2)   = $reads[$i+1] =~ /^.+\/(.+)/;
    my $submit = "";
    if ( $aligner =~ /(?:bowtie2|both)/i )
    {
        my $aln_method    = ( $sample_type =~ /"exome"/i ) ? "--very-sensitive" : "--very-sensitive-local";
        chomp ( my ($btx) = `ls $index_dir/*.4.bt2 | sed s/\.4\.bt2//` );
        submit_step ( "00_aln.bt2",               "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "ALN_METHOD", $aln_method, "BTX", $btx, "R1", $R1, "R2", $R2 );
        submit_step ( "01_sam_to_bam.bt2",        "NAME", $name, "READS_DIR", $reads_dir );
        submit_step ( "02_sort.bt2",              "NAME", $name, "READS_DIR", $reads_dir );
        submit_step ( "03_rm_dups.bt2",           "NAME", $name, "READS_DIR", $reads_dir, "BIN", $bin );
        submit_step ( "04_indel_realigner.bt2",   "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "BIN", $bin, "FASTA", $fasta, "DBSNP", $dbsnp );
        submit_step ( "05_base_recalibrator.bt2", "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "BIN", $bin, "FASTA", $fasta, "DBSNP", $dbsnp );
        submit_step ( "06_genotyper.bt2",         "NAME", $name, "READS_DIR", $reads_dir );
    }
    if ( $aligner =~ /(?:bwa|both)/i )
    {
        submit_step ( "00_aln.bwa",               "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "R1", $R1, "R2", $R2, "FASTA", $fasta );
        submit_step ( "01_sam_to_bam.bwa",        "NAME", $name, "READS_DIR", $reads_dir );
        submit_step ( "02_sort.bwa",              "NAME", $name, "READS_DIR", $reads_dir );
        submit_step ( "03_rm_dups.bwa",           "NAME", $name, "READS_DIR", $reads_dir, "BIN", $bin );
        submit_step ( "04_indel_realigner.bwa",   "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "BIN", $bin, "FASTA", $fasta, "DBSNP", $dbsnp );
        submit_step ( "05_base_recalibrator.bwa", "NAME", $name, "READS_DIR", $reads_dir, "THREADS", $threads, "BIN", $bin, "FASTA", $fasta, "DBSNP", $dbsnp );
        submit_step ( "06_genotyper.bwa",         "NAME", $name, "READS_DIR", $reads_dir );
    }
}

sub usage
{
    die <<USAGE;

    Usage: perl $0 -c <configuration_file.txt>

    Your configuration file MUST have all of these variables defined in it.

    Configuration options available

      VARIABLE    Description                                                   Options
      BIN         Absolute location of the Picard Tools and GATK jar files      Not applicable
      INDEX_DIR   Absolute location of your bowtie2 / bwa indexes               Not applicable
      READS_DIR   Absolute location of the fastq files to be used               Not applicable
      FASTA       Absolute location of the reference fasta file                 Not applicable
      DBSNP       Absolute location of the dbnsp vcf file                       Not applicable
      HAPMAP      Absolute location of the hapmap vcf file                      Not applicable
      OMNI        Absolute location of the omni vcf file                        Not applicable
      MILLS       Absolute location of the mills vcf file                       Not applicable
      ANNOTATOR   The variant annotator to use                                  snpEff, ANNOVAR
      SAMPLE_TYPE The type of sequencing that was done                          exome, rnaseq
      CANCER      Whether or not this is a cancer sample                        yes, no, true, false
      ALIGNER     The aligner software that to use                              bowtie2, bwa, both
      GENOTYPER   The genotyper to use                                          GATK (UnifiedGenotyper), mpileup
      THREADS     Number of threads to use in parallelizable modules            1 < n < max number of cores on the system
      NAME        The name you want to give to this experiment                  Serenity, bob, fancy_pink_ponies, I really don't care; just don't use spaces in the name.

USAGE
}

sub submit_step
{
    my $file = shift;
    my @vars = @_;
    open OUT, ">tmp.pbs";
   #my $script = `cat qsub_files/$file.pbs`;
    my $script = `cat qsub_files/a.pbs`;
    for ( my $i = 0; $i < @vars; $i += 2 )
    { 
        my $a = $vars[$i];
        my $b = $vars[$i+1];
        $script =~ s/$a/$b/g;
    }
    print OUT $script;
   #chomp ( $submit = `qsub -W depend=afterok:$submit tmp.pbs` );
    chomp ( my $submit = `qsub -W depend=afterok:$submit tmp.pbs` );
    close OUT;
}

##
# Assumptions made:
#     index directory has a sym link to the fasta file in that directory or the fasta file is actually in the index directory itself
#         this assumption is made because bwa needs the index file and fasta file to coexist
#     all samples are paired-end
#     the sample name does not have an underscore in it
