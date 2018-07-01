#!/usr/bin/perl
##by Li Lei, 20180629, in St.Paul;
#this is to match the physical position of the SNPs to the Fst values
#usage: ./extraction_SNPs.pl ../data/WBDC_genotype_count_transposed#.txt ../data/all_2H_SNP_namesOnly.txt >all_2H_SNP_WBDC_genotype_count.txt
use strict;
use warnings;

use Data::Dumper;
my ($file1, $file2) = @ARGV;
#fle1: vcf ;file2:Fst

my %gidhash;


open(VCF,  "$file1") or die "Could not open $file1";
my $header = <VCF>;
foreach my $row (<VCF>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $SNPid = $rtemp[2];
        my $Pos = $rtemp[0]."_".$rtemp[1];;
           $gidhash{$SNPid}=$Pos;
        #print "$rtemp[0]\n";
}
close (VCF);

open(SNPID,  "$file2") or die "Could not open $file2";
my $header1 = <SNPID>;
print "SNP\tFST\tChr_2016\tPhysPos_2016\n";
foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $SNP;
        #if ($rtemp[0] =~ /^X*+/){
        	#print "$rtemp[0]\n";
        #    $SNP = $rtemp[0] =~ s/^X//;
         #   print "$SNP\n";
       # }
        #else{
        $SNP = $rtemp[0];
        #	print "$SNP\n";
        #}
        $rtemp[0] =~ s/^X//;
        #print "$rtemp[0]\t$SNP\n";
	    if (exists $gidhash{$rtemp[0]}){
	    	my @phy = split(/\_/,$gidhash{$rtemp[0]});
	    	#print "$SNP\t$rtemp[1]\t$gidhash{$rtemp[0]}\n";
	    	print "$SNP\t$rtemp[1]\t$phy[0]\t$phy[1]\n";
	    }
}
close (SNPID);