#!/usr/bin/perl -w
use strict;

my @a = split /_/, $ARGV[0];
print "\@RG\\tID:$a[0]\\tSM:$a[1]\\tPL:ILLUMINA\\tLB:$a[2]\\tPU:$a[3]";