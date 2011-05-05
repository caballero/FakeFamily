#!/usr/bin/perl

=head1 NAME

parseVCF.pl

=head1 DESCRIPTION

Parse VCF files according to one designated population.

=head1 USAGE

usage: parseVCF.pl POPULATION

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

$ARGV[0] or die "usage: parseVCF.pl POPULATION\n";

my $pop = shift @ARGV;
opendir D, "." or die "cannot open local directory\n";
while (my $file = readdir D) {
    next unless ($file =~ /$pop.vcf.gz/);
    warn "parsing $file\n";
    my $file_h = $file;
    $file_h = "gunzip  -c $file | " if ($file =~ m/gz$/);
    open F, "$file_h" or die "cannot open $file\n";
    while (<F>) {
        next if (m/^#/);
        my @line   = split (/\t/, $_);
        my $chr    = $line[0];
        my $pos    = $line[1];
        my $rs     = $line[2];
        my $ref    = $line[3];
        my $alt    = $line[4];
        my @alt    = split (/,/, $alt);
        my %allele = ();
        my $sum    = 0;
        for (my $i = 9; $i <= $#line; $i++) {
            if ($line[$i] =~ m/(\d+)\/(\d+):\d+/) {
                my $al1 = $1;
                my $al2 = $2;
                
                if ($al1 == 0) { $al1 = $ref; }
                else           { $al1 = $alt[$al1 - 1]; }
            
                if ($al2 == 0) { $al2 = $ref; }
                else           { $al2 = $alt[$al2 - 1]; }
            
                $allele{"$al1/$al2"}++;
                $sum++;
            }
        }
        
        my $allele_freq = undef;
        foreach my $allele (keys %allele) {
            my $freq = $allele{$allele} / $sum;
            $allele_freq .= "$allele=$freq,"; 
        }
        
        $allele_freq =~ s/,$//;
        print join "\t", $chr, $pos, $rs, $ref, $alt, $allele_freq;
        print "\n";
    }
    close F;
}
