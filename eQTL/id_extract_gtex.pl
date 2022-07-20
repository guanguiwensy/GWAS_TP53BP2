#!/usr/bin/perl -w
use strict;
die "perl $0 <list><file>\n" unless  @ARGV==2;
my ($lst,$fi)=@ARGV;
open IN,$lst||die;
my %ha;
map{chomp;$ha{(split)[0]}=1}<IN>;
close IN;

open IN,$fi||die;
my %out;
while(<IN>){
	chomp;
	my $info=(split)[0];
	print $/ if(exists $ha{$info});
}
close IN;


