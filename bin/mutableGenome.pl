#!/usr/bin/perl

=head1 NAME

mutableGenome.pl

=head1 DESCRIPTION

identify regions of a genome that can be mutated.

=head1 USAGE

perl mutableGenome.pl -g genome.fa -e file1.var,file2.var -o output

OPTIONS
    Paramater        Description             Value        Default
    -g --genome      Fasta sequences         FILE*
    -e --exclude     VAR files (to remove)   FILE*
    -o --output      Output file             FILE**       STDOUT
    
    -h --help        Print this screen
    -v --verbose     Verbose mode on
    
     * Files can be compressed with gzip/bzip2, multiple files are
       separated by commas.
    ** Output file is a tad-delimited text:
       CHROM POSITION REF-BASE

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2011

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $output     = undef;
my $genome     = undef;
my $exclude    = undef;
my @exclude_f  = ();
my @genome_f   = ();
my %genome     = ();
my ($genome_h, $exclude_h, $genome_f, $exclude_f);
my ($chr, $pos, $ref, $len);

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'g|genome:s'      => \$genome,
    'e|exclude:s'     => \$exclude,
    'o|output:s'      => \$output
);
pod2usage(-verbose => 2)     if (defined $help);
pod2usage(-verbose => 2) unless (defined $genome);

warn "loading genome\n" if (defined $verbose);
@genome_f = split (/,/, $genome);
foreach $genome_f (@genome_f) {
    $genome_h = $genome_f;
    $genome_h = "gunzip  -c $genome_f | " if ($genome_f =~ m/gz$/);
    $genome_h = "bunzip2 -c $genome_f | " if ($genome_f =~ m/bz2$/);
    open G, "$genome_h" or die "cannot open $genome_f\n";
    while (<G>) {
        chomp;
        if (m/>/) {
            $chr =  $_;
            $chr =~ s/>//;
            warn "    reading $chr\n" if (defined $verbose);
        }
        else {
            $genome{$chr} .= $_;
        }
    }
    close G;
}

warn "masking sites\n" if (defined $verbose);
@exclude_f = split (/,/, $exclude);
foreach $exclude_f (@exclude_f) {
    warn "    using $exclude_f\n" if (defined $verbose);
    $exclude_h = $exclude_f;
    $exclude_h = "gunzip  -c $exclude_f | " if ($exclude_f =~ m/gz$/);
    $exclude_h = "bunzip2 -c $exclude_f | " if ($exclude_f =~ m/bz2$/);
    open E, "$exclude_h" or die "cannot open $exclude_f\n";
    while (<E>) {
        next if (m/#/);
        ($chr, $pos) = split (/\t/, $_);
        next unless ($chr =~ m/[\d+|Y|X]/);
        $chr = "chr$chr";
        substr($genome{$chr}, $pos - 1, 1) = 'X';
    }
    close E;
}

warn "writing mutable genome\n" if (defined $verbose);
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write in $output\n";
}

foreach $chr (sort keys(%genome)) {
    $len = length $genome{$chr};
    warn "    processing $chr ($len bp)\n";
    print ">$chr\n";
    while ($genome{$chr}) {
        my $seq = substr($genome{$chr}, 0, 50);
        print "$seq\n";
        substr($genome{$chr}, 0, 80) = '';
    }
}

