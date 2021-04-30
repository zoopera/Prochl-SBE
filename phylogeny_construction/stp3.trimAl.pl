#!/usr/bin/perl

@methods=("mnf");
for $method (@methods)
{
	system "ls *.$method > list.txt";
	open FILE, "<", "list.txt";
	while(<FILE>)
	{
		chomp;
		if(/^(\S+)\.$method/)
		{
			system "trimal -in $1.$method -out $1.$method.trimAl -gappyout";
		}
	}
	close FILE;
	system "rm list.txt";
}

