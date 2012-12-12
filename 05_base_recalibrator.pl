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
##PBS -N $name.bt2.base_recalibrator
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=32GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bt2.base_recalibrator.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bt2.base_recalibrator.stdout
cd $reads_dir
java -jar $bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $fasta -knownSites $dbsnp -I tmp/$name.indels_realigned.bt2.bam -o tmp/$name.recal_data.bt2.grp
java -jar $bin/GenomeAnalysisTK.jar -T PrintReads -R $fasta -BQSR tmp/$name.recal_data.bt2.grp -I tmp/$name.indels_realigned.bt2.bam -o tmp/$name.BQSR.bt2.bam
perl 06_genotyper.pl -c $config_file -n $name
END
    open OUT, ">qsub/05_$name.bt2.base_recalibrator.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 05_$name.bwa*.qsub" );
}
if ( $aligner =~ /(?:bwa|both)/i )
{
    my $sub_script = <<END;
#!/bin/sh
##PBS -N $name.bwa.base_recalibrator
##PBS -l nodes=1:ppn=$threads
##PBS -l mem=32GB
##PBS -l walltime=24:00:00
##PBS -e /work/unmc_ngs/acornish/$name.bwa.base_recalibrator.stderr
##PBS -o /work/unmc_ngs/acornish/$name.bwa.base_recalibrator.stdout
cd $reads_dir
java -jar $bin/GenomeAnalysisTK.jar -T BaseRecalibrator -R $fasta -knownSites $dbsnp -I tmp/$name.indels_realigned.bwa.bam -o tmp/$name.recal_data.bwa.grp
java -jar $bin/GenomeAnalysisTK.jar -T PrintReads -R $fasta -BQSR tmp/$name.recal_data.bwa.grp -I tmp/$name.indels_realigned.bwa.bam -o tmp/$name.BQSR.bwa.bam
perl 06_genotyper.pl -c $config_file -n $name
END
    open OUT, ">qsub/05_$name.bwa.base_recalibrator.qsub";
    print OUT $sub_script;
    close OUT;
   #system ( "qsub 05_$name.bwa*.qsub" );
}
