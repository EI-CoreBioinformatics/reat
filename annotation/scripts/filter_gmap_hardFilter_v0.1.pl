#!/usr/bin/env perl
#
#	This script is to get the GMAP hits at required identity and coverage cutoff, irrespective of how you run GMAP (with -n0,n1 or n*)
#
#
# AUTHOR: Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk || gemygk@gmail.com)

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $prog = basename($0);
my $usage = "

$prog -- Works for the output of gmap-20130720,gmap-20140804,gmap-20141206,gmap-20141222,gmap-20141229 in gff2_gene format

The input GFF is filtered based on user defined quality and coverage cutoffs. 
The output is to STDOUT in GFF format.

# Any alignment should pass the filter criteria, irrespective of chimeric alignment

Usage: 
	perl $prog --gff <file.gff3>

--gff            -- file.gff3  -- GFF3 file (GMAP output --format gff3_gene) [Required]

--identity       -- [INT 0-100] (>=value Default == 0 )
--coverage       -- [INT 0-100] (>=value Default == 0 )
--keep_original  -- Maintain original GMAP format
                    Otherwise the mRNA line is modified for browser view

--help           -- prints this option menu and quit

\n";

my $gff;
my $cutoff_identity=0;
my $cutoff_coverage=0;
my $keep_original = undef;
my $help;

GetOptions (	"gff=s" => \$gff,
				"identity=i" => \$cutoff_identity,
				"coverage=f" => \$cutoff_coverage,
				"k|keep_original" => \$keep_original,
				"h|help" => \$help,
);

if ($help) {
    die $usage;
}

if (@ARGV) {
    die "\nERROR: Do not understand options: @ARGV\n";
}

unless ($gff) {
    die "\nERROR: No input GFF3 file provided.\n$usage\n";
}

unless(-e $gff) {
		print "\nERROR: Cannot open file $gff:$!\n";
		exit;
	}
my $exon=0;
my %gene_hash=();
my %mrna_hash=();
open(FILE,"< $gff") or die "Cannot open $gff.$!\n";
while (<FILE>) { # Read lines from file(s) specified on command line. Store in $_.
    next if /^#/; # skip comments from $_.
    next unless /\S/; # \S matches non-whitespace.  If not found in $_, skip to next line.
	chomp;
	my @f = split (/\t/); # split $_ at tabs separating fields.
	if($f[2] eq "gene") {
		my ($gene) = $f[8] =~ /ID\s*=\s*([^;]+)/;
		unless (exists $gene_hash{$gene}) {
			$gene_hash{$gene} = $_;
		}
	}
	elsif($f[2] eq "mRNA") {
		my ($mrna) = $f[8] =~ /ID\s*=\s*([^;]+)/;
		my ($parent) = $f[8] =~ /Parent\s*=\s*([^;]+)/;
		my ($curCov) = $f[8] =~ /coverage\s*=\s*([^;]+)/;
		my ($curID) = $f[8] =~ /identity\s*=\s*([^;]+)/;

		my ($name) = $f[8] =~ /Name\s*=\s*([^;]+)/;
		$name =~ s/\s+$//;
		my ($matches) = $f[8] =~ /matches\s*=\s*([^;]+)/;
		$matches =~ s/\s+$//;
		my ($mismatches) = $f[8] =~ /mismatches\s*=\s*([^;]+)/;
		$mismatches =~ s/\s+$//;
		my ($indels) = $f[8] =~ /indels\s*=\s*([^;]+)/;
		$indels =~ s/\s+$//;
		my ($unknowns) = $f[8] =~ /unknowns\s*=\s*([^;]+)/;
		$unknowns =~ s/\s+$//;

		# Main filtering applied here:
		if ($curCov >= $cutoff_coverage && $curID >= $cutoff_identity) {
			# first print gene for pass filter mRNA
			if (exists $gene_hash{$parent}) {
				print $gene_hash{$parent} , "\n"; # No need to delete as this will be just one unique set anyway for GMAP
			}
			# Print the mRNA line and also store the PF mrna to hash
			unless (exists $mrna_hash{$mrna}) {
				$mrna_hash{$mrna} += 1;
			}
			if($keep_original) {
				print "$_\n";
			} else {
				# Print the modified mRNA lines
				print "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[4]\t$f[5]\t$f[6]\t$f[7]\tID=$mrna;Name=$name;Parent=$parent;Note=cov:$curCov|ID:$curID|match:$matches|mismatch:$mismatches|INDEL:$indels|unkwn:$unknowns\n";
			}
		}
	} else {
		my ($parent) = $f[8] =~ /Parent\s*=\s*([^;]+)/;
		# Only print if exists in the mRNA pass filer hash
		if (exists $mrna_hash{$parent}) {
				print "$_\n";
			}
	}
}
exit;