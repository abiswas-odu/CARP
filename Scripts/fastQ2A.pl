#!/usr/bin/perl

use strict;
use warnings;
use File::Spec;
use Bio::Perl;

my ($fileName1, $fileName2) = @ARGV;
my $numSequencesFile1 = numSequences($fileName1);
my $reportEveryFile1 = int($numSequencesFile1/100) || 1;
print "$numSequencesFile1 sequences to convert in $fileName1\n";

my $numSequencesFile2=numSequences($fileName2);
my $reportEveryFile2=int($numSequencesFile2/100) || 1;
print "$numSequencesFile2 sequences to convert in $fileName2\n";

my $in1=Bio::SeqIO->new(-file=>$fileName1,-format=>"fastq-illumina");
my $in1=Bio::SeqIO->new(-file=>$fileName2,-format=>"fastq-illumina");
my $baseSeqFileName = substr($fileName1,0,rindex($fileName1,'.'));
my $seqOut=Bio::SeqIO->new(-file=>">$baseSeqFileName",-format=>"fasta");


while(my $seq=$in1->next_seq){
    $seqOut->write_seq($seq);
}
print "Done with subfile $fileName1.\n";

while(my $seq=$in2->next_seq){
    $seqOut->write_seq($seq);
}
print "Done with subfile $fileName2.\n";


sub numSequences{
  my $file=@_;
  my $num=`wc -l $file`;
  chomp($num);
  $num = $num / 4;
  return $num;
}
