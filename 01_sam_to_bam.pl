#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopt ( "c:n:", \%opt );
my $name         = $opt{n};
my $config       = `cat $opt{c}`;
my $aligner      = ( $config =~ /ALIGNER\s+(\S+)/ )     ? $1 : "bowtie2"; # can be "bwa", "bowtie2", or "both"
my $reads_dir    = ( $config =~ /READS_DIR\s+(\S+)/ )   ? $1 : ".";

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
