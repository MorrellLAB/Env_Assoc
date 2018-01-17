#!/usr/bin/perl
###2018/1/17 by Li Lei
###aim: calculate the minor allele frequency for the 9k datasets Fumi used for GWAS;
##usage: 
use Data::Dumper;

my ($table, $win) = @ARGV;

open (INFILE, "< $table")or die "Can't open $table";

print "SNP\tMaf\n";
my $header = <INFILE>;
while (<INFILE>) {
	my $line = $_;
	chomp $line;
    my($SNP, @tmp) = split(/\t/, $line);
    my $countA=0;
    my $countB=0;
    foreach my $allele (@tmp){
    	if ($allele == 0){
    		$countA++;
    	}
    	elsif($allele == 2){
    		$countB++;
    	}
    }
    my $freA = $countA/($countA+$countB);
    my $freB = $countB/($countA+$countB);
    my $maf = 0;
    if ($freA >= $freB){
        $maf = $freB;
    }
    else{
    	$maf = $freA;
    }
    print "$SNP\t$maf\n";
}
#print Dumper(\%gene_codon);
close(INFILE);
