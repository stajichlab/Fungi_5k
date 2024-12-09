#!/usr/bin/env perl

use strict;
use warnings;

while(<>){
	chomp;
	my @row = split(/,/,$_);
	print ">$row[0]\n",$row[-1],"\n";
}
