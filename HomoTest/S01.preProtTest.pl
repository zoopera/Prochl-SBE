#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

my $dir = "./";
opendir(my $dh,$dir) or die "Can't open $dir: $!";
my @maf = grep {/mnf$/ && -f "$dir/$_"} readdir($dh);
closedir($dh);

foreach my $maf (sort @maf)
{
	my %aln = ();
	my %gap = ();
	my $len = 0;

	my $fasta_in = Bio::SeqIO->new(-file=>"$dir/$maf",-format=>"fasta");
	while(my $seq = $fasta_in->next_seq)
	{
		my ($gnm) = $seq->id =~ /\|(\S+)/;
		$gnm =~ s/_.*//g;
		$len = $seq->length if $len == 0;
		@{$aln{$gnm}} = split(//,$seq->seq);
		for(my $i=0; $i<=$#{$aln{$gnm}}; $i++)
		{
			$gap{$i} = 1 if $aln{$gnm}[$i] eq "-";
		}
	}$fasta_in->close();

	my $trm = $maf;
	$trm =~ s/.*\///g;
	$trm =~ s/maf/rmgap.phy/g;
	my $num = keys %aln;
	open(FH_OUT,">",$trm) or die "Can't open $trm: $!";
	print FH_OUT "$num ",($len - keys %gap),"\n";
	foreach my $gnm (sort keys %aln)
	{
		print FH_OUT $gnm, " " x (15-length($gnm));
		for(my $i=0; $i<=$#{$aln{$gnm}}; $i++)
		{
			print FH_OUT $aln{$gnm}[$i] if !exists $gap{$i};
		}
		print FH_OUT "\n";
	}close(FH_OUT) or die "Can't close $trm: $!";
}


