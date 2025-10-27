#!/usr/bin/env perl

use strict;
use warnings;

while(<>){
	chomp;
	my @row = split(/,/,$_);
	next if $row[0] eq 'intron_id';
	print ">$row[0]\n",$row[-1],"\n";
}
