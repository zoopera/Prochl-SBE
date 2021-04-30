#!/usr/bin/perl

# ===============================================
#  CONCATENAME MSA FROM DIFFERENT PROT FAMILIES
# ===============================================

use strict;
use warnings;

my $dir = ($ARGV[0])? $ARGV[0]:"example";
system "ls $dir/*for_cat > list.txt";

my %seq = ();
my $sid = "";

system("mkdir ptt_fdr.b");

open GENPOS, ">ptt_fdr.b/partition_finder.cfg" or die $!;
# The following parameters should be modified later according to the data
print GENPOS "## ALIGNMENT FILE ##\n";
print GENPOS "alignment = concat.phy;\n\n";
print GENPOS "## BRANCHLENGTHS: linked | unlinked ##\n";
print GENPOS "branchlengths = linked;\n\n";
print GENPOS "## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | <list> ##\n";
print GENPOS "##              for PartitionFinderProtein: all_protein | <list> ##\n";
print GENPOS "models = LG, LG+G, LG+I, LG+I+G, WAG, WAG+G, WAG+I, WAG+I+G, JTT, JTT+G, JTT+I, JTT+I+G;\n\n";
print GENPOS "# MODEL SELECCTION: AIC | AICc | BIC #\n";
print GENPOS "model_selection = BIC;\n\n";
print GENPOS "## DATA BLOCKS: see manual for how to define ##\n";
print GENPOS "[data_blocks]\n";

my @fams = ();

open FILE, "<", "list.txt" or die $!;
while(my $input = <FILE>)
{
	chomp($input);
	
	my $flg = 0;
	my ($fam) = $input =~ /(OG[0-9]+)\./;
	push @fams, $fam;
	print GENPOS "$fam = ";

	open IN, "$input" or die $!;
	while(my $line = <IN>)
	{
		chomp($line);
		if($line =~ /^>(\S+)/)
		{
			$sid = $1;
			$sid =~ s/\|.*//g;
		}else
		{
			my $start = (exists $seq{$sid})? length($seq{$sid})+1:1;
			$seq{$sid} .= $line;
			my $end = length($seq{$sid});
			if($flg == 0)
			{
				print GENPOS "$start-$end;\n";
				$flg = 1;
			}
		}
	}
	close IN;
}
close FILE;

# Continue: partition_finder.cfg
print GENPOS "\n## SCHEMES, search: all | user | greedy ##\n";
print GENPOS "[schemes]\n";
print GENPOS "search = user;\n\n";
print GENPOS "#user schemes go here if search=user. See manual for how to define.#\n";
print GENPOS "together = (",join(", ",@fams),");\n";
print GENPOS "optimal = ";
open IN, "ptt_fdr.a/analysis/best_scheme.txt" or die $!;
while(my $line = <IN>)
{
	chomp($line);
	if($line =~ /Scheme Description in PartitionFinder format/)
	{
		$line = <IN>;
		chomp($line);
		my ($best) = $line =~ /Scheme_step_[0-9]+ = (\(.*\);)/;
		print GENPOS $best,"\n";
		last;
	}
}
print GENPOS "separate = (",join(") (",@fams),");\n";
close GENPOS;

# ==============================================================
# print out the concatinated alignment in fasta format,
# as well as phylip format
# =============================================================

my $lid = 0; # longest sid
my $snm = 0;
my $len = 0;

open CONCAT, ">concat.fas" or die $!;
foreach my $sid (sort{lc($a) cmp lc($b)} keys %seq)
{
	print CONCAT ">$sid\n",$seq{$sid},"\n";
	$lid = length($sid) if length($sid) > $lid;
	$snm++;
	$len = length($seq{$sid});
}
close CONCAT;

open PHYLIP, ">ptt_fdr.b/concat.phy" or die $!;
print PHYLIP "    $snm    $len\n";
foreach my $sid (sort{lc($a) cmp lc($b)} keys %seq)
{
	my $sps = " " x ($lid+5-length($sid));
	print PHYLIP $sid,$sps,$seq{$sid},"\n";
}
close PHYLIP;

system "rm list.txt";
exit 0;



