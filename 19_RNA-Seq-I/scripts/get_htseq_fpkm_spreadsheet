#!/usr/bin/env perl
use warnings;
use strict;


my $USAGE = "usage: perl get_htseq_fpkm_spreadsheet.pl <sample dir> <loc>

where:
<sample dir> is the name of a file with the names of sample directories (no paths)
<loc> is the path to the dir with the sample directories

";

if (@ARGV<2) {
    die $USAGE;
}
unless (-e $ARGV[0]){
    die "\nERROR: file '$ARGV[0]' does not exist.\n\n";
}
unless (-d $ARGV[1]){
    die "\nERROR: directory '$ARGV[1]' does not exist.\n\n";
}

my $LOC = $ARGV[1];
$LOC =~ s/\/$//;
my @fields = split("/", $LOC);
my $last_dir = $fields[@fields-1];
my $study_dir = $LOC;
$study_dir =~ s/$last_dir//;
my $shdir = $study_dir . "shell_scripts";
my $logdir = $study_dir . "logs";
my $htseqdir = $study_dir . "HTSEQ_FPKM";
my $from_dir = "/projects/data/RNASEQ/HTSEQ_FPKM/";
unless (-d $shdir){
    `mkdir $shdir`;
}
unless (-d $logdir){
    `mkdir $logdir`;
}
unless (-d $htseqdir){
    `mkdir $htseqdir`;
}
my $home = $ENV{HOME};
if (-e "$htseqdir/FINAL_master_list_of_genes_counts.htseq.txt"){
    `rm $htseqdir/FINAL_master_list_of_genes_counts.htseq.txt`;
}
my $s = `ln -s $from_dir/FINAL_master_list_of_genes_counts_MIN.GCB535_2016.htseq.txt $htseqdir/FINAL_master_list_of_genes_counts.htseq.txt`;

if (-e "$htseqdir/FINAL_master_list_of_genes_counts.htseq.FPKM.txt"){
    `rm $htseqdir/FINAL_master_list_of_genes_counts.htseq.FPKM.txt`;
}
$s = `ln -s $from_dir/FINAL_master_list_of_genes_counts_MIN.GCB535_2016.htseq.FPKM.txt $htseqdir/FINAL_master_list_of_genes_counts.htseq.FPKM.txt`;
