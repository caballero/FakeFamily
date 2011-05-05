#!/usr/bin/perl

=head1 NAME

parseVCF.pl

=head1 DESCRIPTION

Parse dbSNP VCF file.

=head1 USAGE

usage: parse_dbSNP_VCF.pl dbSNP.vcf

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

$ARGV[0] or die "usage: parse_dbSNP_VCF.pl dbSNP.vcf\n";

my $file = shift @ARGV;
warn "parsing $file\n";
my $file_h = $file;
$file_h = "gunzip  -c $file | " if ($file =~ m/gz$/);
open F, "$file_h" or die "cannot open $file\n";
print "#CHROM\tPOS\tID\tREF\tALT\tGMAF\n";
while (<F>) {
    next if (m/^#/);
    next unless (m/VC=SNP/);
    my @line   = split (/\t/, $_);
    my $chr    = $line[0];
    my $pos    = $line[1];
    my $rs     = $line[2];
    my $ref    = $line[3];
    my $alt    = $line[4];
    next unless (length $ref == 1);
    my $gmaf   = '-';
    if (m/GMAF=(\d\.\d+)/) {
        $gmaf  = $1;
    }
    
    print join "\t", $chr, $pos, $rs, $ref, $alt, $gmaf;
    print "\n";
}
close F;
