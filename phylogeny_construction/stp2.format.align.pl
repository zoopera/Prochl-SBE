#!/usr/bin/perl
$method="mnf";
system "ls *.$method > list.txt";
open FILE, "<", "list.txt";

while (<FILE>)
{
	chomp;
	if (/^(\S+)\.$method$/)
	{
		$fn = $1;
		open IN, "<", "$fn.$method";
		open OUT, ">", "$fn.temp";

		while (<IN>)
		{
			chomp;
			/^\s*$/ and next;
			if (/^>/)
			{
				$_ =~ s/>//;
				print OUT "$seq\n";
				print OUT ">$_\n";
				$seq = '';
			}
			else
			{
				$seq .= $_;
				$seq =~ s/\s//g;
			}
		}
		print OUT "$seq\n";

		close IN;
		close OUT;

		open SEQ, "<", "$fn.temp";
		open SORT, ">", "$fn.sortedseq";

		%hash=();
		while (<SEQ>)
		{
			chomp;
			/^\s*$/ and next;
			if (/^>(.*)$/)
			{
				$seq = <SEQ>;
				chomp $seq;
				$hash{$1} = $seq;
			}
		}

		@sorted = sort {$a cmp $b} (keys %hash);
		foreach $item (@sorted)
		{
			print SORT ">$item\n";
			print SORT "$hash{$item}\n";
		}

		close SEQ;
		close SORT;

		system "rm $fn.temp";
		system "mv $fn.sortedseq $fn.$method";
	}
}

close FILE;
system "rm list.txt";
