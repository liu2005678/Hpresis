#!/usr/bin/perl
use lib "E:/perl_scripts";
use read_fastx;
$file= shift @ARGV;

my (%hash) =read_fastx::read_fasta($file);
$file=~/(Hpfe\d+)/;
$strain=$1;


foreach my $key(keys %hash){
 $hash2{$strain."_".$key}=$hash{$key};
}
$file2=$file."_added".".fa";
$file2=~s/\.fasta//;
open File, ">$file2";
foreach $key(keys %hash2){
	print File ">",$key,"\n",$hash2{$key},"\n";
}