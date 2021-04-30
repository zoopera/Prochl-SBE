#!/usr/bin/env perl

use strict;
use warnings;

# ===================================================================================================
# Get fams showing compositional homogeneity
my $dir = ".";
opendir(my $dh,$dir) or die "Can't open $dir: $!";
my @p4 = grep {/p4$/ && -f "$dir/$_"} readdir($dh);
closedir($dh);

my %fhomo = ();
foreach my $p4 (sort @p4)
{
	my ($fam) = $p4 =~ /(OMA[0-9]+)\.rmgap/;
	open(FH_IN,"<","$dir/$p4") or die "Can't open $p4: $!";
	while(my $line = <FH_IN>)
	{
		chomp($line);
		if($line =~ /All Sequences\s+([0-9\.]+)/)
		{
			$fhomo{$fam} = $1 if $1>=0.05;
		}
	}close(FH_IN) or die "Can't close $p4: $!";
}

# ====================================================================================================
# print summary
open(FH_OUT,">","S05.parP4CmpHomoTest.tbl") or die "Can't open S05.parP4CmpHomoTest.tbl: $!";
foreach my $fam (sort keys %fhomo)
{
	print FH_OUT "$fam\t",$fhomo{$fam},"\n";
}close(FH_OUT) or die "Can't close S05.parP4CmpHomoTest.tbl: $!";
