#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Temp;
use Bio::DB::Sam;


my ($fqfile1, $fqfile2, $out_prefix, $chr_info, $index, $help, $ov_th, $sam, $low_mem);
my $threads = 1;
my $trim = 0;
my $tmp = ".";
#my $flag = "";
#my $bowtie2_path = "/usr/local/bin/bowtie2/bwt2indexes/";
my $bowtie2_path = "/var/data/indexes/bowtie2/";

GetOptions
(
    "f1=s" => \$fqfile1,
    "f2=s" => \$fqfile2,
    "o=s"  => \$out_prefix,
    "s=s"  => \$chr_info,
    "i=s"  => \$index,
    "p=i"  => \$threads,
    "t=i"  => \$trim,
    "sam=s" => \$sam,
    "bp=s" => \$bowtie2_path,
    "ov_th" => \$ov_th,
    "lm"   => \$low_mem,
    "tmp=s"  => \$tmp,
    "h"    => \$help,
);

## Dependencies

my %progs = (
    'wigToBigWig' => 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig',
    'bowtie2' => 'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml',
    'samtools' => 'http://www.htslib.org/',
    'bedSort' => 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedSort',
);

if ($low_mem) {
    $progs{'bedGraphToBigWig'} = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig';
    delete $progs{'wigToBigWig'};
}

foreach(sort keys %progs) {
    unless (`which $_` ne '') {
        die "It seems that $_ is not available. Please install it ($progs{$_}) and place it in PATH.\n";
    }
}

if ($help) {help();}
unless (((($fqfile1 && $fqfile2 && $index) || $sam) && $chr_info && $out_prefix)) {help();}
if ($sam && ($fqfile1 || $fqfile2 || $index)) {
    print "Warning: -sam option selected; -f1, -f2 or -i options will not be used\n";
}
if (!$ov_th && $threads > 6) {
    $threads = 6;
}
if (!$sam && $index) {
        $index  = $bowtie2_path."/".$index;
        unless (check_index($index) == 1) {
            die "Sorry, it seems that $index is not a valid Bowtie2 index. Exiting...\n";
        }
    }

    
my $bamfile = $out_prefix.".bam";
if ($sam) {
    my $temp_log = File::Temp->new();
    system("samtools view -b $sam | samtools sort -@ $threads - | samtools rmdup - $bamfile 2> $temp_log");
}
else {
    print "mapping with Bowtie2...";
    my $bw_log = $out_prefix."_bowtie2.log";
    bowtie($fqfile1, $fqfile2, $index, $bamfile, $threads, $bw_log, $trim, $tmp);
    print " done!\n";
}
print "Converting to BED...";
my $bed = bam_to_bed($bamfile);
print " done!\nGenerating bigWig...";
bed_to_bigWig($bed, $chr_info, $low_mem);
print " done!\nFinished!!";

#########################################################################################
#########################################################################################



sub help {
    print "ATAC_pipe.pl\nThis script will process the ATAC data from a pair of fastq files\n";
    print "to a final bigWig file that can be uploaded to UCSC Browser.\n\nOPTIONS:\n";
    print "\t-f1\t\tfastq file 1. Can be gzipped (required)\n";
    print "\t-f2\t\tfastq file 2. Can be gzipped (required)\n";
    print "\t-o\t\toutput prefix (required)\n";
    print "\t-s\t\tchromosomes info file (required)\n";
    print "\t-i\t\tBowtie2 index, without path (required)\n";
    print "\t-p\t\tnumber of cores to use (default: 1, limited to 6)\n";
    print "\t-t\t\tnumber of bases to be trimmed from each read 3' end (default: 0)\n";
    print "\t-sam\t\talready mapped reads in sam format (directly from bowtie2). Not compatible with -f and -i options\n";
    print "\t-bp\t\tpath to Bowtie2 indexes (default: '/var/data/indexes/bowtie2/')\n";
    print "\t-ov_th\t\toverride the max cores limit (6)\n";
    print "\t-lm\t\tuse less memory when generating bigWig file (will take longer to run)\n";
    print "\t-tmp\t\ttemporal path and prefix for temporal sorting files (default: current dir)\n";
    print "\t-h\t\tshow this message and exit\n\n";
    exit;
}

sub check_index {
    my ($index) = @_;
    $index .= ".1.bt2";
    if (-s $index) {
       return 1;
    }
    else {
        return 0;
    }
}

sub bowtie {
    my ($fq1, $fq2, $index, $output, $cores, $log_file, $trim3, $temp) = @_;
    my $sort_cores = 1;
    if ($cores >= 4) {$sort_cores = 2;}
    my $tmp_log = File::Temp->new();
    system ("echo bowtie2 -t --no-unal --no-mixed -X 2000 -p $cores -3 $trim3 -x $index -1 $fq1 -2 $fq2 > $log_file");
    system("bowtie2 -t --no-unal --no-mixed -X 2000 -p $cores -3 $trim3 -x $index -1 $fq1 -2 $fq2 2>> $log_file |
           samtools view -b - | samtools sort -@ $sort_cores -T $temp - | samtools rmdup - $output 2> $tmp_log");
}

sub bam_to_bed {
    my ($infile) = @_;
    my $sam = Bio::DB::Sam->new(-bam  => $infile);
    my $qual_thres = 10;    # quality threshold
    my $min_length = 130;   # minimum length of fragments (130 for keeping only nucleosome free)
    my $bedfile = file_name($infile)."_nucfree.bed";
    my @flags = (99, 147, 83, 163);
    open(BED, ">", $bedfile) or die "Could not open $bedfile for writing. Aborting...\n";

    my $iter = $sam->get_seq_stream();
    while (my $read = $iter->next_seq) {
        my $len = $read->isize;
        unless(($read->qual >= $qual_thres) && (abs($len) <= $min_length) && (grep { $_ eq $read->flag } @flags)) {next};
        my $chr = $read->seq_id;
        my ($start, $strand);
        if ($read->strand == 1) {
            $start = $read->start - 1;  ## +4 - 5
	    $strand = '+';
        }
        elsif ($read->strand == -1) {
            $start = $read->end - 10;  ## -5 - 5
	    $strand = '-';
        }
        else {next;}
        my $end = $start + 10;
        
        print BED "$chr\t$start\t$end\t$len\t0\t$strand\n";
     }

    return $bedfile;
}

sub bed_to_bigWig {
    my ($infile, $sizes, $lm) = @_;
    my $file_out = file_name($infile).".bw";
    my $track_file = file_name($infile)."_track.txt";
    my $track_name = file_name($infile);
    if ($lm) {
        my $tmp_bg = File::Temp->new();
        system("bedSort $infile stdout | genomeCoverageBed -bg -i - -g $sizes > $tmp_bg");
        system("bedGraphToBigWig $tmp_bg $sizes $file_out");
    }
    else{
        system("bedSort $infile stdout | genomeCoverageBed -bg -i - -g $sizes | wigToBigWig stdin $sizes $file_out");
    }
    open(TRACK, ">", $track_file) or die "Could not open $track_file: $!\n";
    print TRACK "track type=bigWig name='$track_name' description='$track_name' color=0,0,0 visibility=2 maxHeightPixels=128:60:11 bigDataUrl=http://193.147.188.155/tracks/ATAC/$file_out\n";
    close TRACK;
}

sub file_name {
    my ($file) = @_;
    my @fname = split /\./, $file;
    pop @fname;
    return join ".", @fname;
}
