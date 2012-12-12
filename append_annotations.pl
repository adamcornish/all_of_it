#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Std;

my %opt;
my @db   = "";
my $curr = "";
getopts ( "i:", \%opt );

print "Usage: perl $0 -i all.vcf.\n\tSimply provide the vcf file that contains the combined annotated results.\n\n";
my $file = $opt{i};
# read in the vcf file
my @data = `grep -vP '^#' $file`;
# get the names of all of the samples
chomp ( my $tmp = `grep -P '^#CHROM' $file` );
$tmp =~ s/#CHROM.+?FORMAT\s//;
my @names = split /\t/, $tmp;
open OUT, ">all.variants.txt";
# regex
my $important = "(?:NON_SYNONYMOUS_CODING|NON_SYNONYMOUS_START|STOP_GAINED|STOP_LOST|SYNONYMOUS_CODING|SYNONYMOUS_STOP)";
# set up the initial table
print OUT "Chromosoms\tPosition\tRef\tAlt\tQuality\tTotal Coverage\tGene Symbol\tExon ID\tEffect\tImpact\tOld_codon/New_codon\tOld_AA/New_AA\tdbSNP ID\tPhyloP\tSIFT\tPolyphen2\tLRT\tMutationTaster\tCOSMIC\t";
for ( my $i = 0; $i < @names; $i++ )
{
    my $name = $names[$i];
    print OUT "Ref depth-$name\tAlt depth-$name\t% depth-$name";
    ( $i == $#names) ? print OUT "\n" : print OUT "\t";
}

# Start looping through the lines in the vcf file
foreach my $datum ( @data )
{
    my @split = split /\t/, $datum;
    $split[2] =~ s/\./N\\a/;
    my $meta = $split[7];
    my ($eff) = $split[7] =~ /EFF=(\S+)/;
    my @effects = split /,/, $eff;
    my ($DP) = $datum =~ /DP=(\d+)/;
    my $pre = "$split[0]\t$split[1]\t$split[3]\t$split[4]\t$split[5]\t$DP\t";
    my $post = "";
    my $chr = $split[0];
    if ( $curr ne $chr )
    {
        print "\tReading in $chr.\n";
        @db = `cat /safer/genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbNSFP/$chr | grep -vP '^#'`;
        $curr = $chr;
    }
    for ( my $i = 9; $i < @split; $i++ )
    {
        if ( $split[$i] =~ /\.\/\./ )
        {
            $post .= "0\t0\t0.00%";
        }
        else
        {
            my ( $ref, $alt) = $split[$i] =~ /.+?:(\d+),(\d+)/;
            my $percent = "0\%"; # assume you will divide by 0
            $percent = sprintf ( "%.2f%%", ($alt/($ref+$alt))*100 ) unless ( $ref + $alt == 0 );
            $post .= "$ref\t$alt\t$percent";
        }
        $post .= ( $i == $#split) ? "\n" : "\t";
    }
    foreach my $effect ( @effects )
    {
        if ( $effect =~ /$important/ )
        {
            my ($type, $impact, $codon_change, $AA_change, $gene_symbol, $exon_ID) = $effect =~ /(.+?)\((.*?)\|.*?\|(.*?)\|(.*?)\|.*?\|(.*?)\|.*?\|.*?\|.*?\|(.*?)\)/;
            print OUT "$pre$gene_symbol\t$exon_ID\t$type\t$impact\t$codon_change\t$AA_change\t$split[2]\t$post";
            my ($loc, $parens) = $effect =~ /(.+?)\((.+?)\)/;
            my ($impact,$type,$codon,$aa,$gene,$rna,$coding,$acc,$exon) = split /\|/, $parens;
            print OUT "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$loc\t$impact\t$type\t$codon\t$aa\t$gene\t$acc\t$exon";
            if ( $effect =~ /NON_SYN/)
            {
                my ($ref_AA, $alt_AA) = $aa =~ /([A-Z])\d+([A-Z])/;
                my ($ref_nuc, $alt_nuc, $index) = ($split[3], $split[4], $split[1]);
                my $result = bin_search ( \@split, \@db, $ref_AA, $alt_AA );
                if ( $result )
                {
                    my @col = split "\t", $result;
                    my $out = "\t$col[7]\t$col[8]\t$col[9]\t$col[10]\t$col[12]";
                    print OUT $out;
                }
            }
        }
        print OUT "\n";
    }

    close OUT;
}

sub bin_search
{
    my ( $arr, $db, $ref_AA, $alt_AA ) = @_;
    my $max = $#$db;
    my $min = 0;

    while ( $max >= $min ) 
    {
        my $mid = int( ( $max + $min ) / 2 );
        my @s   = split "\t", $$db[$mid];
        if    ( $s[6] < $$arr[1] ) { $min = $mid + 1; }
        elsif ( $s[6] > $$arr[1] ) { $max = $mid - 1; }
        else
        { 
            my $regex = "$alt_AA\\s$$arr[1]";
            for ( my $i = ($mid - 4); $i < ($mid + 4); $i++ )
            {
                return $$db[$i] if $$db[$i] =~ /$regex/;
                return "" if ($i == ($mid + 3));
            }
        }
    }
}
