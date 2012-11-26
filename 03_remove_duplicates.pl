#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

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
##PBS -N $name.bt2.mark_duplicates
##PBS -l select=1
##PBS -l mem=8GB
##PBS -l walltime=8:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.mark_duplicates.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.mark_duplicates.stdout
cd $reads_dir
java -jar $bin/MarkDuplicates.jar I=tmp/$name.bt2.sorted.bam O=tmp/$name.bt2.dup_removed.bam REMOVE_DUPLICATES=true M=tmp/$name.bt2.mark_dups_metrics_file
perl 04_indel_realigner.pl -c $config_file -n $name
END
    open OUT, ">qsub/03_$name.bt2.mark_duplicates.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 03_$name.bt2*.qsub" );
}

if ( $aligner =~ /(?:bwa|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa.mark_duplicates
##PBS -l select=1
##PBS -l mem=8GB
##PBS -l walltime=8:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.mark_duplicates.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.mark_duplicates.stdout
cd $reads_dir
java -jar $bin/MarkDuplicates.jar I=tmp/$name.bwa.sorted.bam O=tmp/$name.bwa.dup_removed.bam REMOVE_DUPLICATES=true M=tmp/$name.bwa.mark_dups_metrics_file
perl 04_indel_realigner.pl -c $config_file -n $name
END
    open OUT, ">qsub/03_$name.bwa.mark_duplicates.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 03_$name.bwa*.qsub" );
}
