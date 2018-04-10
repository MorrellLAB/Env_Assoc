#!/usr/bin/perl
###2018/4/3 by Li Lei
###aim: match the ancestral SNP with the genotyping file and get the dereived sfs;
##usage: 
use Data::Dumper;

my ($nam, $table) = @ARGV;

open (INFILE, "< $nam")or die "Can't open $nam";

my $header = <INFILE>;
my %hash;
while (<INFILE>) {
    my $line = $_;
    chomp $line;
    my @tmp = split(/\t/, $line);
       $hash{$tmp[0]} = $tmp[1];
}
#print Dumper(\%gene_codon);
close(INFILE);



open (INFILE, "< $table")or die "Can't open $table";

my $header = <INFILE>;
print "SNP\tCHR\tPOS\tREF\tALT\tDAF\n";
while (<INFILE>) {
    my $line = $_;
    chomp $line;
    my($SNP, $chr, $pos, $ref, $alt, @tmp) = split(/\t/, $line);
    my $id = $chr."_".$pos;
    #count the reference and alternative alleles;
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
    my $derive_f;

    if (exists $hash{$id}){
        my $a_a = $hash{$id};
        if ($ref eq $a_a){
            $derive_f = $freB;
        }
        elsif ($alt eq $a_a){
            $derive_f = $freA;
        }
        else{
            $derive_f = "NA";
        }    
    }
    else{
        $derive_f = "NA";
    }

    print "$SNP\t$chr\t$pos\t$ref\t$alt\t$derive_f\n";
    
}
#print Dumper(\%gene_codon);
close(INFILE);