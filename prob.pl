#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use File::Basename; #library for filename separation
# use POSIX;

# require 'dumpvar.pl'; # see dumpval
$| = 1; # turn on autoflush

################################################################################
use constant USAGE =><<END;
Usage: $0 [OPTIONS]
    Reads a sequence file with an RNA sequence; supported formats are 
        FASTA, Vienna, GenBank, and EMBL. 
    Determines pairing probabilities at a given temperature via
        RNAfold and RNAplfold.
    Extracts from RNA*fold's output the pairing probabilities, and 
    plot it via GLE.

    OPTIONS:
    -f <gene.vie[.gz]>     Input file (lines starting with [!%#] are not read)
                                T is equivalent to U; [^AUCG]->N
   [-m "regex"]            cleavage site; allows for sequence or regexp
                                (only containing square bracket pairs []); 
                                default: "GU"
   [-o <#,#,...>]          Start, end of ORF(s)
   [-t <#>]                Temperature in degree celsius (20 <= t <= 50)
                                (default: 37 degree Celsius; uses best 
                                thermodynamic parameters)
   [-p <gene.prob>]        Output file with pairing probabilities
   [-g <gene_RNAfold.gle>] Output file (GLE script for fold data)
   [-u <#>]                Mean prob of unpaired regions (l=1 to l=u)
                               for RNAplfold (default 10)
   [-W <#>]                Window size for RNAplfold (default 80)
   [-L <#>]                Separation of a base pair for RNAplfold (default 40)
   [-k]                    Keep the intermediate files of RNA*fold (default: no)
   [-H <header>]           Printed above the graphics (ONE word!) 

END
################################################################################
my $cmdline = basename($0); 

my %opts;
getopts('hvkm:f:t:p:g:o:u:W:L:H:', \%opts);

# help
if ($opts{h}) { print USAGE; exit;}

# verbose; only for debugging
my $verbose = 0;    # no
if ($opts{v}) {
    $verbose = 1;
}
my $debug = 0;  # only for the strange things in RNAplfold

my $Header = "";
if ($opts{H}) {
    $Header = $opts{H};
}

# keep intermediate files of RNA*fold?
my $keep = 0;   # no
if ($opts{k}) {
    $keep = 1;
}

# sequence file
my $infile = "";
if (!$opts{f}) {
    die USAGE . "\nERROR: No file with RNA sequence given.\n";
} else {
    $infile = $opts{f};
    if (!-e $infile) {
        die USAGE . "\nERROR: $infile does not exist.\n";
    }
}
# remember file name parts to create output file names
my ($in_id, $in_dir, $in_ext) = fileparse($infile, '\..*');
#if ($verbose) {
#    print "in_ID  = $in_id\n";
#    print "in_DIR = $in_dir\n";
#    print "in_EXT = $in_ext\n";
#}

# Output file with pairing probabilities
my $probfile = "";
if (!$opts{p}) {
    $probfile = $in_dir . $in_id . ".prob";
} else {
    $probfile = $opts{p};
}

# GLE script and data file
my $glefile = "";
my $datfile = "";
if (!$opts{g}) {
    $glefile = $in_dir . $in_id . ".gle";
    $datfile = $in_dir . $in_id . ".dat";
} else {
    $glefile = $opts{g};
    my ($id, $dir, $ext) = fileparse($glefile, '\..*');
    $datfile = $dir . $id . ".dat";
}

# temperature; 37 degree celsius uses best thermodynamic parameters
my $temp = "37";
if ($opts{t}) {
    $temp = $opts{t};
    if ($temp<20 || $temp>50) {
        die USAGE . "\nERROR: Temperature outside range 20 <= t <= 50.\n";
    }
}

# comma-separated list of sequence positions (start1,end1,start2,end2,...) 
#   to label, f. e., ORFs;
#   all are plotted in yellow; that is, different ORFs can't be differentiated 
#       if regions overlap
my @orf;
if ($opts{o}) {
    @orf = split(/,/, $opts{o});
}

# Window size for RNAplfold
my $winSize = 80;
if ($opts{W}) {
    $winSize = $opts{W};
        if ($winSize<20 || $winSize>500) {
        die USAGE . "\nERROR: Window size outside range 20 <= W <= 500.\n";
    }
}

# Span size (max separation of a basepair) for RNAplfold
my $spanSize = 40;
if ($opts{L}) {
    $spanSize = $opts{L};
        if ($spanSize<20 || $spanSize>500) {
        die USAGE . "\nERROR: Span size outside range 20 <= L <= 500.\n";
    }
}

# Mean prob of unpaired regions for RNAplfold
my $upregion  = 10;
my $upregion2 =  5;
if ($opts{u}) {
    $upregion = $opts{u};
    if ($upregion<5 || $upregion>50) {
        die USAGE . "\nERROR: Length of unpaired region outside range 5 <= u <= 50.\n";
    }
    $upregion2 = $upregion/2; # slight ERROR if upregion isn't even
}

# cleavage site; allows for sequence or simple regexp
my $motif = "GU";
if ($opts{m}) {
    $motif = uc($opts{m});
}

# log used options in GLE's datfile
$cmdline .= 
    " -f " . $infile   . 
    " -m " . $motif;
if ($opts{o}) {
    $cmdline .= " -o " . $opts{o};
}
$cmdline .= 
    " -t " . $temp     .
    " -p " . $probfile .
    " -g " . $glefile  .
    " -u " . $upregion .
    " -W " . $winSize  .
    " -L " . $spanSize .
    " -H " . $Header;
    
if ($verbose) {
    print "infile   = $infile\n";
    print "probfile = $probfile\n";
    print "glefile  = $glefile\n";
    print "datfile  = $datfile\n";
    print "motif    = $motif\n";
    print "temp     = $temp °C\n";
    print "W size   = $winSize\n";
    print "L size   = $spanSize\n";
    print "u size   = $upregion\n";
    if (exists $orf[0]) {
        print "ORF(s) = ";
        for (my $i=0; $i<@orf; $i+=2) {
            if (!exists $orf[$i+1]) {
                die "You have to give an even number of ORF limits!\n";
            }
            print $orf[$i]."--".$orf[$i+1].', '
        }
        print "\n";
    } else {
        print "No ORFs given\n";
    }
}

################################################################################
# Check for executables
use File::Which;
my $pgm_RNAfold   = which "RNAfold";
my $pgm_RNAplfold = which "RNAplfold";
# my $pgm_gle       = which "gle";

if ($verbose) {
    print "RNAfold:   $pgm_RNAfold\n";
    print "RNAplfold: $pgm_RNAplfold\n";
#     print "gle:       $pgm_gle\n";
}

my $y_min = 1e-2;   # min of prob axis

################################################################################
# return max value of array
sub max {
    my ($max) = shift(@_);
    my $foo;
    foreach $foo (@_) {
        $max = $foo if $max < $foo;
    }
    return $max;
}
# return min value of array
sub min {
    my ($min) = shift(@_);
    my $foo;
    foreach $foo (@_) {
        $min = $foo if $min > $foo;
    }
    return $min;
}

################################################################################
# Determine sequence file format
# Return either FASTA(=Vienna), GenBank or EMBL
sub seqFileFormat {
    my $file = shift;

	my $format = "";
    open(INP, "$file") or die "Can't open \"$file\": $!";
	    while (<INP>) {
	        $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	        chomp;
	        # print "FORMAT: $_\n";
            if      (m/^>\s*.*$/) {
                $format = "FASTA";
            } elsif (m/^LOCUS\s+/) {
                $format = "GB";
            } elsif (m/^ID\s+/) {
                $format = "EMBL";
            } else {
                die USAGE . "\nERROR: Unknown file format!\n";
            }
            last;
        }
    close(INP);
	# print "FORMAT: $format\n";
    return $format;
}
################################################################################
# Read sequence from file
# 
sub readFASTA {
	my $file = shift;

    my $seq;
	my $header;
	open(INP, "$file") or die "Can't open \"$file\": $!";
	    while (<INP>) {
	        $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	        chomp;
	        if (m/>\s*(.*)$/) {
	            $header = $1;
	            $seq    = "";
	        } elsif (m/^([ACGTNURYMWKSBDHV]+)$/i) {
	            $seq .= $1;
	        } elsif (m/\s*/) {
	            # empty line
	        } elsif (m/^[#!%]/) {
	            # comment line
	        } else {
	            die "ERROR FASTA: Strange line in input $infile: >$_<\n";
	        }
	    }
	close(INP);
	$seq = uc($seq);
	$seq =~ tr/TRYMWKSBDHV/UNNNNNNNNNN/;
    return ($header, $seq);
}
sub readGB {
	my $file = shift;

    my $seq = "";
	my $header = "";
    my $last = 0;
	open(INP, "$file") or die "Can't open \"$file\": $!";
	    while (<INP>) {
	        if ($last==1) {last;}
            $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	        chomp;
	        if (m/^LOCUS\s+(\w+)/) {
	            # print "GB HEADER: $1\n";
	            $header = $1;
	        } elsif (m/^ORIGIN\s*$/) {
                while (<INP>) {
	                $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	                chomp;
	                # print "GB seq: $_\n";
                    if (m/\s*\d+\s+([ACGTNURYMWKSBDHV ]+)$/i) {
                        $seq .= $1;
                    } elsif (m,^//$,) {
	                    # print "GB SEQ: $seq\n";
                        $last = 1;
                        last;
                    }
                }
	        }
	    }
	close(INP);
	$seq = uc($seq);
	$seq =~ tr/TRYMWKSBDHV/UNNNNNNNNNN/;
    $seq =~ s/\s//g;
    return ($header, $seq);
}
sub readEMBL {
	my $file = shift;

    my $seq = "";
	my $header = "";
    my $last = 0;
	open(INP, "$file") or die "Can't open \"$file\": $!";
	    while (<INP>) {
	        if ($last==1) {last;}
	        $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	        chomp;
	        if (m/^ID\s+(\w+)/) {
	            $header = $1;
	        } elsif (m/^SQ\s+/) {
                while (<INP>) {
	                $_ =~ s/\r\n?/\n/;  # Windoof: LF/CR; Linux: only CR
	                chomp;
	                if (m/\s+([ACGTNURYMWKSBDHV ]+)\s+\d+$/i) {
                        $seq .= $1;
                    } elsif (m,^//$,) {
                        $last = 1;
                        last;
                    }
                }
	        }
	    }
	close(INP);
	$seq = uc($seq);
	$seq =~ tr/TRYMWKSBDHV/UNNNNNNNNNN/;
    $seq =~ s/\s//g;
    return ($header, $seq);
}

$infile =~ s/(.*\.gz)\s*$/gzip --decompress --to-stdout $1 |/;
my $format = seqFileFormat($infile);
my $header = "";
my $seq    = "";
if      ($format eq "FASTA") {
    ($header, $seq) = readFASTA($infile);
} elsif ($format eq "GB") {
    ($header, $seq) = readGB($infile);
} elsif ($format eq "EMBL") {
    ($header, $seq) = readEMBL($infile);
}
my $length = length($seq);
if ($length==0) {
    die USAGE . "\nERROR: No sequence read\n";
}

if ($verbose) {
    print "HEADER:>$header<\n";
    print "SEQ:>$seq<\n";
    print "SEQ len: $length\n";
}

my $dummy = $header;
   $dummy = substr($dummy,0,29);    # RNAfold uses first 29 chars for filenames
   $dummy =~ s/ .*//;               #   but truncates at a blank
my $rnafold_dp_file = $dummy . "_dp.ps";
my $rnafold_ss_file = $dummy . "_ss.ps";
my $rnafold_lunp    = $dummy . "_lunp";
if ($verbose) {
    print "DP file  = $rnafold_dp_file\n";
    print "SS file  = $rnafold_ss_file\n\n";
}
# write sequence in Vienna format
my $tmpfile = $header . ".vie";
open(VIE,">$tmpfile") or die "ERROR: Can't write sequence file!\n";
    print VIE ">$header\n$seq\n";
close(VIE);

### determine positions of GUC in sequence
#   store pos of G(uc) in @pos_GUC
#   store range (-10/+12) around G(uc) in @range_GUC
my @pos_GUC    = ("*") x ($length+1);
my @range_GUC  = (0)   x ($length+1);
my $copy_motif = $motif;
my $motif_length = 0;
# print "motif        = $motif\n";
if ($copy_motif =~ m/\[/) { # contains regex; do not allow for more complex things then []
    $motif_length  = ( $copy_motif =~ s/\[.*?\]//g);    # count bracket pairs
    $motif_length += length($copy_motif);               # add remaining nucleotids if any
} else {
    $motif_length  = length($copy_motif);               # simple string
}
# print "motif_length = $motif_length\n";
while ($seq =~ m/($motif)/g) {
#my $bla = pos($seq);
#print "position     = $bla\n";
#print "motif_length = $motif_length\n";
#my $blabla = $bla-$motif_length+1;
#print "pos          = " . $blabla . "\n";
#print "seq          = " . substr($seq, $blabla, 5) . "\n";
#exit;
    $pos_GUC[pos($seq)-$motif_length+1] = $y_min*2;
    for (my $i=max(1,pos($seq)-$motif_length+1-10); $i<=min(pos($seq)+10,$length); $i++) {
        $range_GUC[$i] = 1;
    }
}

################################################################################
# RNAfold
# 
my $fold_options = '-d2 --noLP ';
if ($temp != 37) {
    $fold_options .= '-T ' . $temp . ' ';
}

my $rnafold = $pgm_RNAfold . ' ' . $fold_options . '-p  < ' . $tmpfile . ' |';
if ($verbose) {
    print "$rnafold\n";
}
my $logF = "";
open (RNAFOLD, "$rnafold") or die "Can't open \"$rnafold\": $!";
    while (<RNAFOLD>) {
        chomp;
        if (m/^scaling factor/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^free energy/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^[\.*]+/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^ frequency/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^>/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^[ACGUNacgun]+$/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^\s*$/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } elsif (m/^\s*[()\{\}\|,.]+/) {
            if ($verbose) {print $_ . "\n";}
            $logF .= "!" . $_ . "\n";
        } else {
            print "Strange output from RNAfold:\n\t$_\nProceed with care!\n";
        }
    }
close(RNAFOLD);

################################################################################
# read output of RNAfold
# 
my $rnafold_seq = "";
my @mm = (0) x ($length+1); # nts paired in mfe structure
my @pp = (0) x ($length+1); # pairing prob from partition function
my @sp = (0) x ($length+1); # sp[i] = -Sum p_i * ln(p_i) = measure of well-definedness
my $meanP = 0;
open(RNAFOLD, "<$rnafold_dp_file") or die "Can't open \"$rnafold_dp_file\": $!";
    while (<RNAFOLD>) {
        chomp;
        if (/\/sequence \{ \((\S*)[\\\)]/) {
        	$rnafold_seq = $1;              # empty for new version
	        while (!/\) \} def/) {  # read until end of definition
	            $_ = <RNAFOLD>;
	            /(\S*)[\\\)]/;      # ends either with `)' or `\'
	            $rnafold_seq .= $1;
	        }
	        next;
        }

        next unless /(\d+) (\d+) (\d+\.\d+) (.box)$/;
        my ($i, $j, $p, $id) = ($1,$2,$3,$4);
        if ($id eq "ubox") {
	        $p *= $p;           # square it to probability
            $meanP +=  $p;
	        my $ss = $p>0 ? $p*log($p) : 0;
	        $sp[$i] += $ss;
	        $sp[$j] += $ss;
	        $pp[$i] += $p;
	        $pp[$j] += $p;
        }
        if ($id eq "lbox") {
	        $mm[$i]++;
	        $mm[$j]++;
        }
    }
close(RNAFOLD);
$meanP /= $length;

my $log2 = log(2);
for (my $i=1; $i<=$length; $i++) {
    if ($pp[$i]<1.0) {
        $sp[$i]  = -1*((1-$pp[$i])*log(1-$pp[$i]))/$log2;
    } elsif ($pp[$i]>1.0) {
        printf "STRANGE RNAfold: p(%d) = %5.3f\n", $i, $pp[$i];
    }
}

if (uc($rnafold_seq) ne uc($seq)) {
    print "Input seq and RNAfold seq differ:\n$seq\n$rnafold_seq\nProceed with care!\n\n";
}

################################################################################
# RNAplfold
# 
my $rnaplfold = $pgm_RNAplfold . ' ' . $fold_options . ' -W ' . $winSize . ' -L ' . $spanSize . ' -u ' . $upregion .  ' < ' . $tmpfile . ' |';
if ($verbose) {
    print "$rnaplfold\n";
}
my $logL = "";
open (RNAPLFOLD, "$rnaplfold") or die "Can't open \"$rnaplfold\": $!";
    while (<RNAPLFOLD>) {
        chomp;
        if (m/^>/) {
            if ($verbose) {print $_ . "\n";}
            $logL .= "!" . $_ . "\n";
        } else {
            print "Strange output from RNAplfold:\n\t$_\nProceed with care!\n";
        }
    }
close(RNAPLFOLD);

################################################################################
# read output of RNAplfold
# 
my @ppL = (0) x ($length+1); # pairing prob from partition function
my @spL = (0) x ($length+1); # sp[i] = -Sum p_i * ln(p_i) = measure of well-definedness
my $meanPL = 0;
open(RNAPLFOLD, "<$rnafold_dp_file") or die "Can't open \"$rnafold_dp_file\": $!";
    while (<RNAPLFOLD>) {
        chomp;
        if (/\/sequence \{ \((\S*)[\\\)]/) {
        	$rnafold_seq = $1;      # empty for new version
	        while (!/\) \} def/) {  # read until end of definition
	            $_ = <RNAPLFOLD>;
	            /(\S*)[\\\)]/;      # ends either with `)' or `\'
	            $rnafold_seq .= $1;
	        }
	        next;
        }

        if (m/(\d+) (\d+) (\d+\.\d+) (.box)$/) {
            my ($i, $j, $p, $id) = ($1,$2,$3,$4);
            if ($id eq "ubox") {
	            $p *= $p;           # square it to probability
              # $p  = $p/$length*$winSize;  # Why are there prob values > 1 ?????
                $meanPL += $p;
	            my $ss = $p>0 ? $p*log($p) : 0;
	            $spL[$i] += $ss;
	            $spL[$j] += $ss;
	            $ppL[$i] += $p;
	            $ppL[$j] += $p;
            }
        } elsif (m/^\/winSize (\d+) def/) {
            if ($spanSize != $1) {die "ERROR in RNAplfold: win/spanSize $spanSize != $1\n";}
        }
   }
close(RNAPLFOLD);
$meanPL /= $length;

for (my $i=1; $i<=$length; $i++) {
    if ($ppL[$i]<1.0) {
        $spL[$i]  = -1*((1-$ppL[$i])*log(1-$ppL[$i]))/$log2;
    } elsif ($ppL[$i]>1.0) {
        if ($debug) {printf "STRANGE RNAplfold: p(%d) = %5.3f\n", $i, $ppL[$i];}
    }
}

my @lunp; # stores prob that [i-$upregion+1..i] is unpaired
# print "rnafold_lunp = $rnafold_lunp\n";
open(RNAPLFOLD, "<$rnafold_lunp") or die "Can't open \"$rnafold_lunp\": $!";
    while (<RNAPLFOLD>) {
        chomp;
        if (m/^(\d+)\s+(.*)/) {
            my $index = $1;
            my $dummy = (split(/\s+/,$2))[-1];
            $dummy =~ s/NA/*/g;
            $lunp[$index] = $dummy;
        }
    }
close(RNAPLFOLD);
# for (my $i=1; $i<scalar(@lunp); $i++) {
#     printf "%4d\t", $i;
#     for (my $j=0; $j<scalar(@{$lunp[$i]}); $j++) {
#         print $lunp[$i][$j] . "\t";
#     }
#     print "\n";
# }
################################################################################

open(GLEDAT, ">$datfile") or die "Can't open \"$datfile\": $!";
    print  GLEDAT '! ' . $cmdline . "\n";
    my @nt = split("", $seq);
    print  GLEDAT '!' . $rnafold . "\n";
    print  GLEDAT       $logF;
    print  GLEDAT '!' . $rnaplfold . "\n";
    print  GLEDAT       $logL;
    print  GLEDAT "!  pp[i] = pairing prob from partition function\n";
    print  GLEDAT "!  mm[i] = nts paired in mfe structure\n";
    print  GLEDAT "!  sp[i] = -Sum pp[i] * ln(pp[i]) = measure of well-definedness\n";
    print  GLEDAT "! pfL[i] = pairing prob from local partition function\n";
    print  GLEDAT "! spL[i] = -Sum ppL[i] * ln(ppL[i]) = measure of well-definedness\n";
    printf GLEDAT "! %s    = position of %s in sequence\n", $motif, $motif;
    printf GLEDAT "! %s<   = position of %s in sequence with pp[i]<%f\n", $motif, $motif, $y_min*10;
    printf GLEDAT "! %s<<  = position of %s in sequence with pp[i]<%f\n", $motif, $motif, $y_min;
    printf GLEDAT "! lunp   = mean probability that regions of length 1 to length %d centered at i are unpaired\n", $upregion;
    my $i_fmt = length(sprintf("%d", $length));
    printf GLEDAT "!%2s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%5s\t%s\t%s\t%s\t%s\t%s\n", 
                    "i", "pp","pp","mm","sp","pfL","pfL","spL", $motif, $motif."<".$y_min*10, $motif."<".$y_min, "lunp", "seq";
    for (my $i=1; $i<=$length; $i++) {
        printf GLEDAT "%${i_fmt}d\t", $i;
        if ($range_GUC[$i] == 0) {
            printf GLEDAT "%-5s\t", "1.0";
        } else {
            printf GLEDAT "%5.3f\t", $pp[$i];
        }
        printf GLEDAT "%5.3f\t%5.3f\t%5.3f\t", $pp[$i], $mm[$i], $sp[$i];
        if ($range_GUC[$i] == 0) {
            printf GLEDAT "%-5s\t", "1.0";
        } else {
            printf GLEDAT "%5.3f\t", $ppL[$i];
        }
        printf GLEDAT "%5.3f\t%5.3f\t", $ppL[$i], $spL[$i];
        if ($pos_GUC[$i] eq "*") {
            print  GLEDAT " *   \t";
        } else {
            printf GLEDAT "%5.3f\t", $pos_GUC[$i];
        }

#        if ($pos_GUC[$i] ne "*" && $pp[$i]<$y_min*10) {
#            printf GLEDAT "%f\t", $pos_GUC[$i];
#        } else {
#            print  GLEDAT "*\t";
#        }
#        if ($pos_GUC[$i] ne "*" && $pp[$i]<$y_min) {
#            printf GLEDAT "%f\t", $pos_GUC[$i];
#        } else {
#            print  GLEDAT "*\t";
#        }
#
#   1 2 3 4 5 6 7 8 9 10
#   +-+-+-+-+-+-+-+-+-+
#           G U C
#  i-4      i
#        printf "%4d\t%s\t", $i, $pos_GUC[$i];
#        printf "%s\t", $lunp[$i];
#        printf "%s\t", $lunp[$i+1];
#        printf "%s\t", $lunp[$i+2];
#        printf "%s\t", $lunp[$i+3];
#        printf "%s\t", $lunp[$i+4];
#        printf "%s\n", $lunp[$i+5];
        if ($i+$upregion2<$length && $lunp[$i+$upregion2] ne "*") {
            if (($pos_GUC[$i] ne "*") && ($lunp[$i+$upregion2]>$y_min)) {
                printf GLEDAT "%5.3f\t", $pos_GUC[$i];
            } else {
                print  GLEDAT " *   \t";
            }
        } else {
                print  GLEDAT " *   \t";
        }
        if ($i+$upregion2<$length && $lunp[$i+$upregion2] ne "*") {
            if (($pos_GUC[$i] ne "*") && ($lunp[$i+$upregion2]>$y_min*10)) {
                printf GLEDAT "%5.3f\t", $pos_GUC[$i];
            } else {
                print  GLEDAT " *   \t";
            } 
        } else {
                print  GLEDAT " *   \t";
        }

        if ($range_GUC[$i] == 0) {
            print  GLEDAT " *   \t";
        } else {
            if ($i+$upregion2<$length && $lunp[$i+$upregion2] ne "*") {
                if ($lunp[$i+$upregion2] eq "*") {
                    printf GLEDAT "%-5s\t", "*";
                } else {
                    printf GLEDAT "%5.3f\t", $lunp[$i+$upregion2];
                }
            } else {
                    printf GLEDAT "%-5s\t", "*";
            }
        }
        printf GLEDAT "%s", $nt[$i-1];
        print  GLEDAT "\n";
    }
close(GLEDAT);

################################################################################
################################################################################
# constants used in  GLE file
use constant GLESTART =><<GSTART;
size 29.7 7.0
finalsize = 18/0.7
factor = pagewidth()/finalsize
sub val2cm val unit\$
   if      unit\$ = "inch"   then
      return val*2.54*factor
   else if unit\$ = "bp"     then
      return val/72*2.54*factor
   else if unit\$ = "pt"     then
      return val/72.27*2.54*factor
   else if unit\$ = "pica"   then
      return val/6.0225*2.54*factor
   else if unit\$ = "didot"  then
      return val/67.553*2.54*factor
   else if unit\$ = "cicero" then
      return val/5.6294*2.54*factor
   else if unit\$ = "cm" then
      return val*factor
   end if
end sub
titlehei = val2cm(12,"bp")
charhei  = val2cm(10,"bp")
smallhei = val2cm( 8,"bp")
foothei  = val2cm( 7,"bp")
tinyhei  = val2cm( 6,"bp")
set hei charhei
lsswidth = val2cm(0.25,"bp")
 lswidth = val2cm(0.50,"bp")
 lnwidth = val2cm(1.00,"bp")
 ldwidth = val2cm(2.00,"bp")
lddwidth = val2cm(4.00,"bp")
orange\$   = "rgb255(238,127,0)"
skyeblue\$ = "rgb255(181,203,214)"
green\$    = "rgb255(140,177,16)"
yellow\$   = "rgb255(240,228,66)"
blue\$     = "rgb255(0,106,179)"
red\$      = "rgb255(190,10,38)"
purple\$   = "rgb255(204,121,167)"
set font pshb
sub labelC x y height var\$
   gsave
      set hei height 
      amove xg(x) yg(y)
      rmove -.5*twidth(var\$) -.5*theight(var\$)
      write var\$
   grestore
end sub
sub labelR x y height var\$
   gsave
      set hei height 
      amove xg(x) yg(y)
      rmove -1.*twidth(var\$) -.5*theight(var\$)
      write var\$
   grestore
end sub
sub labelL x y height var\$
   gsave
      set hei height
      amove xg(x) yg(y)
      rmove 0 -.5*theight(var\$)
      write var\$
   grestore
end sub
GSTART
use constant GLEMAIN =><<'GMAIN';
set hei charhei
amove -2 0
begin graph
    title Header$ hei titlehei
    ytitle "Prob"               hei titlehei
    xaxis  min x_min  max x_max hei charhei
    yaxis  min y_min  max y_max hei charhei
    xaxis  dticks delta_x  dsubticks delta_x/2
!    yaxis  dticks delta_y  dsubticks delta_y/2
    yaxis log
    xside                     lwidth lnwidth
    xticks    length lddwidth lwidth lnwidth
    xsubticks length ldwidth  lwidth lswidth
    yside                     lwidth lnwidth
    yticks    length lddwidth lwidth lnwidth
    ysubticks length ldwidth  lwidth lswidth
!    x2axis off
!    y2axis off
GMAIN
################################################################################
use constant GLEEND =><<"GEND";
begin key
   hei foothei
   position tr offset -.1 0
   text key1\$ lstyle 1 lwidth ldwidth color red\$
   text key5\$ lstyle 1 lwidth ldwidth color blue\$
   text key8\$  marker star msize smallhei color grey
   text key9\$  marker star msize charhei
   text key10\$ marker star msize titlehei color green\$
   text \" \"
   text key11\$ line lwidth lswidth color black
   text key0\$
end key
GEND
################################################################################
sub texclean {
    my $string = shift;
    $string =~ s/_/\\_/g;
    return $string;
}
################################################################################
# write GLE file
open(GLE, ">$glefile") or die "Can't open \"$glefile\": $!";
    print GLE GLESTART;
    printf GLE 
"x_min   =  0
x_max   =  %d
delta_x =  %d
y_min   =  %f
y_max   =  1.0
Header\$ = \"%s\"\n",  ceil($length/(10**($i_fmt-2)))*(10**($i_fmt-2)),
                   floor($length/(10**($i_fmt-1)))*(10**($i_fmt-2)), 
                   $y_min,
                   texclean($Header);
    print GLE GLEMAIN;
    printf GLE 
"    xtitle \"Sequence " . texclean($header) . "\"   hei titlehei
    data \"$datfile\" ignore 11
    d1  line lwidth lswidth color red\$         ! i, p (filtered +/-10 pos around %s)
!   d2 line lwidth lswidth color red\$          ! i, p
!   d3 line lwidth lswidth color grey30         ! i, m (mfe)
!   d4 line lwidth lswidth color grey10         ! i, sp (well-defindness of p)
    d5  line lwidth lswidth color blue\$        ! i, p (filtered +/-10 pos around %s) from PLfold
!   d6  line lwidth lswidth color blue\$        ! i, p                                 from PLfold
!   d7  line lwidth lswidth color grey10        ! i, sp (well-defindness of p)         from PLfold
    d8  marker star msize smallhei color grey   ! %s
    d9  marker star msize charhei               ! %s && p<%1.0e
    d10 marker star msize titlehei color green\$! %s && p<%1.0e
    d11 line lwidth lnwidth color black         ! u[i]...u[i+%d]
end graph
key1\$  = \"\\setfont{pshbo}p\\setfont{pshb}_{global}\"
key5\$  = \"\\setfont{pshbo}p\\setfont{pshb}_{local}\"
key8\$  = \"%s\"
key9\$  = \"%s with \\setfont{pshbo}u\\setfont{pshb}(%d) > %1.0e\"
key10\$ = \"%s with \\setfont{pshbo}u\\setfont{pshb}(%d) > %1.0e\"
key11\$ = \"\\setfont{pshbo}u\\setfont{pshb}(%d)\"
key0\$  = \"\\setfont{pshbo}W\\setfont{pshb} = %d; \\setfont{pshbo}L\\setfont{pshb} = %d\"
\n", $motif, 
     $motif, 
     $motif, 
     $motif, $y_min*10, 
     $motif, $y_min, 
     $upregion, 
     $motif, 
     $motif, $upregion, $y_min, 
     $motif, $upregion, $y_min*10, 
     $upregion, $winSize, $spanSize;
    print GLE GLEEND;
    if (exists $orf[0]) {
        print GLE "set lwidth lddwidth\nset color yellow\$\n";
        for (my $i=0; $i<@orf; $i+=2) {
            printf GLE "amove xg(%d) yg(ygmin)+yg(y_min*1.5)-yg(y_min)\n", $orf[$i];
            printf GLE "aline xg(%d) yg(ygmin)+yg(y_min*1.5)-yg(y_min)\n", $orf[$i+1];
        }
    }
    printf GLE "
set color red\$
set lwidth lswidth
amove xg(xgmin) yg(%f)
aline xg(xgmax) yg(%f)
set color blue\$
amove xg(xgmin) yg(%f)
aline xg(xgmax) yg(%f)
set color black
amove xg(%d) yg(ygmin)
aline xg(%d) yg(ygmax)
\n", $meanP, $meanP, $meanPL, $meanPL, $length, $length;
close(GLE);

if (!$keep) {
    unlink($rnafold_dp_file) or die "Can't delete $rnafold_dp_file: $!\n";
    unlink($rnafold_ss_file) or die "Can't delete $rnafold_ss_file: $!\n";
    unlink($rnafold_lunp)    or die "Can't delete $rnafold_lunp: $!\n";
    unlink($tmpfile)         or die "Can't delete $rtmpfile: $!\n";
}
