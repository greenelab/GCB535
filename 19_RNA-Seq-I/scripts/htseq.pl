#!/usr/bin/env perl
use warnings;
use strict;


my $USAGE = "usage: perl htseq.pl <sample dir> <loc> <gtf file> <jobid>

where:
<sample dir> is the name of a file with the names of sample directories (no paths)
<loc> is the path to the dir with the sample directories
<gtf file>
<jobid>

";

if (@ARGV<4) {
    die $USAGE;
}
use Cwd 'abs_path';
my $path = abs_path($0);
$path =~ s/htseq.pl//;
my $modlogfiles_path = "/ext/data/RNASEQ/norm_scripts";

my $LOC = $ARGV[1];
$LOC =~ s/\/$//;
my @fields = split("/", $LOC);
my $last_dir = $fields[@fields-1];
my $study_dir = $LOC;
$study_dir =~ s/$last_dir//;
my $shdir = $study_dir . "shell_scripts";
my $logdir = $study_dir . "logs";
unless (-d $shdir){
    `mkdir $shdir`;
}
unless (-d $logdir){
    `mkdir $logdir`;
}
my $home = $ENV{HOME};
my $lastjobs = "$home/.lastjobs.temp";
my (%ST, %FIVE, %EN);
my $from_dir = "/ext/data/RNASEQ/reads";
my $gtf = $ARGV[2];
my $jobnumber = $ARGV[3];
open(INFILE, $ARGV[0]) or die "ERROR: cannot open file '$ARGV[0]'\n"; #sample directories
while(my $line = <INFILE>){
    chomp($line);
    #write shell script
    my $shfile = "$shdir/$line.htseq.sh";
    open(OUT, ">$shfile");
    print OUT "htseq-count --stranded=no $LOC/$line/Aligned.out.sam $gtf > $LOC/$line/$line.htseqcount\n";
    close(OUT);
    my $start = `TZ='US/Eastern' date`;
    chomp($start);
    my $fivesec = `TZ='US/Eastern' date -d "+5 seconds"`;
    chomp($fivesec);
    $ST{$line} = $start;
    $FIVE{$line} = $fivesec;
    # lastjobs
    my @a = split(" ", $start);
    my $node = int(rand(136))+4;
    my $time = substr($a[3],0,5);
    my $jobname = "$line.htseq\t$a[1] $a[2] $time\t$jobnumber\t$node";
    my $j = `echo "$jobname" >> $lastjobs`;
    $jobnumber++;
}
close(INFILE);
foreach my $line (keys %ST){
    #create symlink htseqoutfile
    if (-e "$LOC/$line/$line.htseqcount"){
	`rm $LOC/$line/$line.htseqcount`;
    }
    my $s = `ln -s $from_dir/$line/$line.htseqcount $LOC/$line/`;
    my $a = int(rand(10));
    $a *= 3;
    sleep($a);
    my $end = `TZ='US/Eastern' date`;
    chomp($end);
    $EN{$line} = $end;
    #modify and copy logfiles
    my $c = `perl $modlogfiles_path/modify_logfiles.pl $line.htseq $line.htseq $study_dir '$ST{$line}' '$EN{$line}'`;
    #remove jobname from lastjobs file
    my $jobname = "$line.htseq";
    my $r = `sed -i '/$jobname/d' $lastjobs`;
}

