#!/usr/bin/env perl
use warnings;
use strict;

my $USAGE = 

"usage: runall_htseq.pl <sample dir> <loc> <gtf>

where:
<sample dir> is the name of a file with the names of sample directories (no paths)
<loc> is the path to the dir with the sample directories
<gtf> full path to the gtf file

-h : display usage

";

if (@ARGV<3){
    die $USAGE;
}

for(my $i=0; $i<@ARGV; $i++){
    if ($ARGV[$i] eq '-h'){
	die $USAGE;
    }
}
use Cwd 'abs_path';
my $path = abs_path($0);
$path =~ s/runall_htseq.pl//;
my $script = "$path/htseq.pl";
my $home = $ENV{HOME};
my $lastjobs = "$home/.lastjobs.temp";

my $LOC = $ARGV[1];
unless (-d "$LOC"){
    die "ERROR: Directory '$LOC' does not exist.\n";
}
my $gtf = $ARGV[2];
unless (-e "$gtf"){
    die "ERROR: gtf file '$gtf' does not exist.\n";
}

open(IN, $ARGV[0]);
while(my $line = <IN>){
    chomp($line);
    unless (-d "$LOC/$line/"){
	die "\nERROR: Directory '$LOC/$line' does not exist.\n\n";
    }
}
close(IN);

if (-e $lastjobs){
    open(IN, $ARGV[0]) or die "ERROR: cannot open file '$ARGV[0]'\n";
    while(my $line = <IN>){
	chomp($line);
	my $name_to_check = "$line.htseq";
	if (`grep -c $name_to_check $lastjobs` > 1){
	    die "ERROR: htseq job for sample $line is already running.";
	}
    }
    close(IN);
}
my $pid = open my $fhOut, "| sleep 1", or die;
chomp($pid);
system("perl $script $ARGV[0] $LOC $gtf $pid&");
