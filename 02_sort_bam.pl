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
my ($bin)       = $config =~ /BIN\s+(\S+)/;
my ($aligner)   = $config =~ /ALIGNER\s+(\S+)/; # can be "bwa", "bowtie2", or "both"
my ($reads_dir) = $config =~ /READS_DIR\s+(\S+)/;

if ( $aligner =~ /(?:bowtie2|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bt2.sort_bam
##PBS -l select=1
##PBS -l mem=4GB
##PBS -l walltime=6:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.sort_bam.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.sort_bam.stdout
cd $reads_dir
java -jar $bin/SortSam.jar INPUT=tmp/$name.bt2.bam OUTPUT=tmp/$name.bt2.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
perl 03_remove_duplicates.pl -c $config_file -n $name
END
    open OUT, ">qsub/02_$name.bt2.sort_bam.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 02_$name.bt2*.qsub" );
}
if ( $aligner =~ /(?:bwa|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa.sort_bam
##PBS -l select=1
##PBS -l mem=4GB
##PBS -l walltime=6:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.sam_to_bam.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.sam_to_bam.stdout
cd $reads_dir
java -jar $bin/SortSam.jar INPUT=tmp/$name.bwa.bam OUTPUT=tmp/$name.bwa.sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
perl 03_remove_duplicates.pl -c $config_file -n $name
END
    open OUT, ">qsub/02_$name.bwa.sort_bam.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 02_$name.bwa*.qsub" );
}
