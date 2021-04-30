#!/usr/bin/perl

$seqFile = $ARGV[0];
$datFile = $ARGV[1];
$fam = $ARGV[2];

# get number of genomes
open IN2, "<", $seqFile;
@array=<IN2>;
$num=(scalar @array - 1)/2; # number of orthologs in gene family $1
close IN2 ;

open IN, "<", $datFile;
open OUT, ">", "$datFile.fm";
while(<IN>)
{
	chomp;
	# format the header line of "Proportions ... (below diagonal) item1 item2 ..."
	if(/^(Proportions.*\(below diagonal\))\s+(\S.*) $/)
	{
		print $datFile;
		print OUT "$1\n";
		@data=split /\s+/, $2;
		$count=0;
		for $item (@data)
		{
			$count++;
			if($count==1)
			{
				print OUT "   $item";
			}
			elsif($count==$num)
			{
				print OUT " \n";
				$count=0;
			}
			else
			{
				print OUT "    $item";
			}
		}
	}
	# otherwise print lines
	else
	{
		print OUT "$_\n";
	}
}
close IN;
close OUT;

system "mv $datFile.fm $datFile";

