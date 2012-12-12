#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
getopt ( "c:n:", \%opt );
my $name        = $opt{n};
my $config      = `cat $opt{c}`;
my ($bin)       = $config =~ /BIN\s+(\S+)/;
my ($aligner)   = $config =~ /ALIGNER\s+(\S+)/;     # can be "bwa", "bowtie2", or "both"
my ($fasta)     = $config =~ /FASTA\s+(\S+)/;
my ($threads)   = $config =~ /THREADS\s+(\S+)/;
my ($reads_dir) = $config =~ /READS_DIR\s+(\S+)/;
my ($dbsnp)     = $config =~ /DBSNP\s+(\S+)/;

if ( $aligner =~ /(?:bowtie2|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bt2.indel_realigner
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=32GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.indel_realigner.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.indel_realigner.stdout
cd $reads_dir
java -jar $bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $fasta -I tmp/$name.bt2.dup_removed.bam -known $dbsnp -o tmp/$name.indel.bt2.intervals
java -jar $bin/GenomeAnalysisTK.jar -T IndelRealigner         -R $fasta -I tmp/$name.bt2.dup_removed.bam -known $dbsnp -o tmp/$name.indels_realigned.bt2.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals tmp/$name.indel.bt2.intervals
perl 05_base_recalibrator.pl -c $config_file -n $name
END
    open OUT, ">qsub/04_$name.bt2.indel_realigner.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 04_$name.bt2*.qsub" );
}
if ( $aligner =~ /(?:bwa|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa.indel_realigner
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=32GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.indel_realigner.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.indel_realigner.stdout
cd $reads_dir
java -jar $bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $fasta -I tmp/$name.bwa.dup_removed.bam -known $dbsnp -o tmp/$name.indel.bwa.intervals
java -jar $bin/GenomeAnalysisTK.jar -T IndelRealigner         -R $fasta -I tmp/$name.bt2.dup_removed.bam -known $dbsnp -o tmp/$name.indels_realigned.bt2.bam --maxReadsForRealignment 100000 --maxReadsInMemory 1000000 -targetIntervals tmp/$name.indel.bt2.intervals
perl 05_base_recalibrator.pl -c $config_file -n $name
END
    open OUT, ">qsub/04_$name.bwa.indel_realigner.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 04_$name.bwa*.qsub" );
}
