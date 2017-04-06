#!/usr/bin/perl
##by Li Lei, 20170406, in St.Paul;
#this is to extract the physical positions for the SNPs Fuimi identified as significant association with the enviromental factors!
#usage: perl script/find_new_physical_pos.pl 9k_SNPs/sorted_all_9k_masked_90idt.vcf sig_ass_SNP/sig_asso_SNPs_fumi.txt  >sig_ass_SNP/sig_asso_SNPs_physical_pos_fumi.txt
use strict;
use warnings;
use Data::Dumper;
my ($file1, $file2) = @ARGV;

my %SNPposHash;


open(SNPID,  "$file1") or die "Could not open $file1";

foreach my $row (<SNPID>){
       if ($row =~ /^#/) {
           next;
       }
       else{
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $SNP_id = $rtemp[2];
           $SNPposHash{$SNP_id} = $rtemp[0]."_".$rtemp[1];
      }
        #print "$rtemp[0]\n";
}
close (SNPID);
#print Dumper(\%gidhash);

open(SIGSNP,  "$file2") or die "Could not open $file2";
my $header = <SIGSNP>;
print "Phenotype\tSNP\tChromosome\tPosition\tP.value\tmaf\tnobs\tRsquare.of.Model.without.SNP\tRsquare.of.Model.with.SNP\tFDR_Adjusted_P.values\n";
foreach my $row (<SIGSNP>){
        chomp $row;
        my @rtemp1 = split(/\t/,$row);
           #print "$rtemp1[0]\n";
           my $SNP = $rtemp1[1];
        if(exists $SNPposHash{$SNP}){
           my $value = $SNPposHash{$SNP};
           my @temp = split(/\_/,$value);
           print "$rtemp1[0]\t$rtemp1[1]\t$temp[0]\t$temp[1]\t$rtemp1[4]\t$rtemp1[5]\t$rtemp1[6]\t$rtemp1[7]\t$rtemp1[8]\t\n";
        }
        else{
           print "$SNP\tNOT_EXISTS\n";
        }
}
close (SNPID);