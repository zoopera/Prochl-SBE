#!/usr/bin/perl

$seqFile = $ARGV[0];
$datFile = $ARGV[1];
$outFile = $ARGV[2];

# get number of genomes
open IN2, "<", $seqFile;
@array=<IN2>;
$num=(scalar @array - 1)/2; # number of orthologs in gene family $1
close IN2 ;

open IN, "<", $datFile;
open OUT, ">>", $outFile;
while(<IN>)
{
	chomp;
	# format the header line of "Proportions ... (below diagonal) item1 item2 ..."
	if(/^(Proportions.*\(below diagonal\))\s+(\S.*) $/)
	{
		print OUT "$datFile\n";
	}
}
close IN;
close OUT;
