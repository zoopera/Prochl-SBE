#!/usr/bin/perl
system "ls *.mnf.trimAl > list.txt";
open FILE, "<", "list.txt";
open TAXA, "<", "genome.list";
open OUTSUM, ">", "summary.missing.gene";

@taxa;
%taxa_count;
while(<TAXA>)
{
	chomp;
	if(/^(\S+)$/)
	{
		push @taxa, $1;
		$taxa_count{$1}=0;
	}
}

while(<FILE>)
{
	chomp;
	if(/^(\S+)\.mnf\.trimAl$/)
	{
		$fam=$1;
		open IN, "<", "$1.mnf.trimAl";
		open OUT, ">", "$1.for_cat";
		%hash=();
		$len=0;
		while(<IN>)
		{
			chomp;
			if(/^>(\S+)\|\S+$/)
			{
				$id=$1;
				$seq=<IN>;
				$len=length $seq;
				chomp $seq;
				$hash{$id}=$seq;
			}
		}
		for $taxon (@taxa)
		{
			if(exists $hash{$taxon})
			{
				$taxa_count{$taxon}++;
				print OUT ">$taxon\n$hash{$taxon}\n";
			}
			else
			{
				$gap='';
				for($x=1; $x<$len; $x++)
				{
					$gap .= '-';
				}
				print OUT ">$taxon\n$gap\n";
			}
		}
		close IN;
		close OUT;
	}
}

for $key (keys %taxa_count)
{
	print OUTSUM "$key\t$taxa_count{$key}\n";
}

close FILE;
close TAXA;
close OUTSUM;
system "rm list.txt";
system "sort -k2n summary.missing.gene > summary.missing.gene2";
system "mv summary.missing.gene2 summary.missing.gene";


