#!/usr/bin/env perl

use strict;
use warnings;

# input: (a) xx.phy (b) xx.tab (prottest)
# output: (a) xx.nex (b) xx.tre (c) xx.py (p4 script)

my $dir = ".";
opendir(my $dh,$dir) or die "Can't open $dir: $!";
my @phy = grep {/phy$/ && -f "$dir/$_"} readdir($dh);
closedir($dh);

foreach my $phy (sort @phy)
{
	my ($fam) = $phy =~ /(\S+).phy/;

	# (a) xx.nex
	my $nex = "$fam.nex";
	my %seq = ();
	open(PHY_IN,"<","$dir/$phy") or die "Can't open $phy: $!";
	while(my $line = <PHY_IN>)
	{
		chomp($line);
		if($line =~ /^([\S]+)\s+([A-Z\-]+)/)
		{
			my $gnm = substr($1,0,10);
			$seq{$gnm} = $2;
		}
	}close(PHY_IN) or die "Can't close $phy: $!";

	my @gnm = sort keys %seq;
	open(NEX_OUT,">",$nex) or die "Can't open $nex: $!";
	print NEX_OUT "#NEXUS\n\n";
	print NEX_OUT "begin taxa;\n";
	print NEX_OUT "  dimensions ntax=",scalar(@gnm),";\n";
	print NEX_OUT "  taxlabels ",join(" ",@gnm),";\nend;\n\n";
	print NEX_OUT "begin characters;\n";
	print NEX_OUT "  dimensions nChar=",length($seq{$gnm[0]}),";\n";
	print NEX_OUT "  format datatype=protein gap=- missing=?;\n";
	print NEX_OUT "  matrix\n";
	foreach my $gnm (@gnm)
	{
		print NEX_OUT "    $gnm";
		my @base = split(//,$seq{$gnm});
		for(my $i=0; $i<=$#base; $i++)
		{
			if($i%60 == 0)
			{
				print NEX_OUT "\n"," " x 6;
			}
			print NEX_OUT $base[$i];
		}
		print NEX_OUT "\n";
	}
	print NEX_OUT "  ;\nend;\n";
	close(NEX_OUT) or die "Can't close $nex: $!";

	# (b) xx.tre
	# (c) xx.py (p4)
	my $ready4parm = 0;
	my $ready4tree = 0;
	my $prtt = "$fam.tab";
	my $tree = "$fam.tre";
	my $py = "$fam.py";
	open(PRTT_IN,"<",$prtt) or die "Can't open $prtt: $!";

	open(TREE_OUT,">",$tree) or die "Can't open $tree: $!";
	print TREE_OUT "#NEXUS\n\n";
	print TREE_OUT "begin taxa;\n";
	print TREE_OUT " dimensions ntax=",scalar(@gnm),";\n";
	print TREE_OUT " taxlabels ",join(" ",@gnm),";\nend;\n\n";
	print TREE_OUT "begin trees;\ntree random = [&U] ";

	open(PY_OUT,">",$py) or die "Can't open $py: $!";
	print PY_OUT "var.doCheckForDuplicateSequences = 0\n";
	print PY_OUT "read('$nex')\n";
	print PY_OUT "a = var.alignments[0]\n";
	print PY_OUT "read('$tree')\n";
	print PY_OUT "t = var.trees[0]\n";
	print PY_OUT "d = Data()\n";
	print PY_OUT "t.data = d\n";
	print PY_OUT "t.newComp(free=1,spec='empirical')\n";

	my $suffix = "";
	while(my $line = <PRTT_IN>)
	{
		chomp($line);
		if($line =~ /Best model according to BIC:\s+(\S+)/)
		{
			my $model = $1;
			$suffix = $1;
			my @suffix = split(/\+/,$suffix);
			$suffix = "+".join("+",@suffix[1..$#suffix]) if scalar(@suffix)>1;
			$suffix = "" if scalar(@suffix)==1;
			if($model =~ /Dayhoff/){$model="d78"};
			if($model =~ /JTT/){$model = "jtt"};
			if($model =~ /LG/){$model = "lg"};
			if($model =~ /RtREV/){$model = "rtRev"};
			if($model =~ /WAG/){$model = "wag"};
			print PY_OUT "t.newRMatrix(free=0,spec='$model')\n";
			print PY_OUT "t.setNGammaCat(nGammaCat=4)\n" if $suffix =~ /\+G/;
			$suffix =~ s/\+F//g;
			$suffix =~ s/\+/\\+/g;
		}elsif($ready4parm == 1 && $line =~ /alpha\s+\($suffix\):\s+([0-9\.]+)/)
		{
			print PY_OUT "t.newGdasrv(free=1,val=$1)\n";
		}elsif($ready4parm == 1 && $line =~ /p-inv\s+\($suffix\):\s+([0-9\.]+)/)
		{
			print PY_OUT "t.setPInvar(free=1,val=$1)\n\n";
		}elsif($line =~ /Sets included in the consensus tree/)
		{
			$ready4tree = 1;
		}elsif($line =~ /Model-averaged estimate of parameters/)
		{
			$ready4parm = 1;
		}elsif($ready4tree == 1 && $line =~ /^\(/)
		{
			print TREE_OUT "$line\nend;\n";
		}elsif($line =~ /Best model according to LnL:/)
		{
			$ready4parm = 0;
			$ready4tree = 0;
			last;
		}
	}
	print PY_OUT "t.setPInvar(free=1,val=0)\n\n" if $suffix !~ /\+I/;
	print PY_OUT "t.compoTestUsingSimulations(nSims=1000,doChiSquare=True)\n";
	close(PRTT_IN) or die "Can't close $prtt: $!";
	close(TREE_OUT) or die "Can't close $tree: $!";
	close(PY_OUT) or die "Can't close $py: $!";
}
