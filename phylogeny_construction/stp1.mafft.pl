#!/usr/bin/perl
system "ls *.pepseq > list.mnf";
open FILE, "<", "list.mnf";

while (<FILE>)
{
	chomp;
	if (/^(\S+)\.pepseq$/)
	{
			open SCRIPT, ">", "$1.mnf.in";
			printf SCRIPT "einsi $1.pepseq > $1.mnf";
			close SCRIPT;
			system "chmod 755 $1.mnf.in";
open SCRIPT, ">", "$1.lsf";
print SCRIPT "#!/bin/sh\n";
print SCRIPT "APP_NAME=AMD_small\n";
print SCRIPT "NP=1\n";
print SCRIPT "RUN=\"RAW\"\n";
print SCRIPT "./$1.mnf.in\n";
close SCRIPT;
system "chmod 755 $1.lsf";
system "bsub $1.lsf";


	}
}

close FILE;
system "rm list.mnf";
