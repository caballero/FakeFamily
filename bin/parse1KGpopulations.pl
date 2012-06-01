#!/usr/bin/perl

=head1 NAME

parse1KGpopulations.pl

=head1 DESCRIPTION

Fetch variation calls (VCFs) from the 1000genomes project, parse the file,
compute allele frequencies per population.

=head1 USAGE

parse1KGpopulations.pl SAMPLE_POPS GENOTYPES.vcf.gz  

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

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

$ARGV[2] or die "use parse1KGpopulations.pl SAMPLE_POPS GENOTYPES.vcf.gz OUT\n";
my $sampop_file = shift @ARGV;
my $vcf_file    = shift @ARGV;
my $out         = shift @ARGV;
my %res;
my %sam2pop;
my %poptot;
my %allpop;
my @allpop;
my @label;
my @line;
my ($pop, $var, $sam, $chr, $pos, $rs, $ref, $alt);
my $min_share = 3;

warn "loading population per sample\n";
open POP, "$sampop_file" or die "cannot open $sampop_file\n";
while (<POP>) {
    chomp;
    ($pop, $sam)   = split (/\t/, $_);
    $sam2pop{$sam} = $pop;
    $allpop{$pop}  = 1;
}
close POP;
@allpop = sort keys %allpop;

warn "parsing VCF file\n";
open VCF, "gunzip -c $vcf_file | " or die "cannot open $vcf_file\n";
while (<VCF>) {
    chomp;
    my %var = ();
    if (m/^#CHROM/) {
        @label = split (/\t/, $_);
        for (my $i = 9; $i <= $#label; $i++) {
            $poptot{ $sam2pop{$label[$i]} }++;
        }
    }
    elsif (m/^#/) {
        next; # skip the rest of the header
    }
    else {
        @line = split (/\t/, $_);
        $chr  = $line[0];
        $pos  = $line[1];
        $rs   = $line[2];
        $ref  = $line[3];
        $alt  = $line[4];
        
        for (my $i = 9; $i <= $#line; $i++) {
            $pop        = 'NA';
            $pop        = $sam2pop{ $label[$i] } if (defined $sam2pop{ $label[$i] });
            ($var)      = split (/:/, $line[$i]);
            $var{$pop}{$var}++;
        }
        foreach $pop (@allpop) {
            next unless (defined $poptot{$pop});
            next if ($poptot{$pop} < $min_share);
            my $f1      = 0; # freq for 0|0 - homo ref
            my $f2      = 0; # freq for 0|1 - hetero
            my $f3      = 0; # freq for 1|0 - hetero
            my $f4      = 0; # freq for 1|1 - homo noref
            if (defined $var{$pop}{'0|0'}) {
                next if ($var{$pop}{'0|0'} == $poptot{$pop}); # all is homo ref
                $f1     = sprintf("%.8f", $var{$pop}{'0|0'} / $poptot{$pop});
            }
            
            my $noref   = 0;
             
            if (defined $var{$pop}{'0|1'}) {
                $noref += $var{$pop}{'0|1'};
                $f2     = sprintf("%.8f", $var{$pop}{'0|1'} / $poptot{$pop});
            }
             
            if (defined $var{$pop}{'1|0'}) {
                $noref += $var{$pop}{'1|0'};
                $f3     = sprintf("%.8f", $var{$pop}{'1|0'} / $poptot{$pop});
            }
             
            if (defined $var{$pop}{'1|1'}) {
                $noref += $var{$pop}{'1|1'};
                $f4     = sprintf("%.8f", $var{$pop}{'1|1'} / $poptot{$pop});
            }
            next if ($noref < $min_share);
            
            $res{"$out.$pop.data"} .= join "\t", $chr, $pos, $rs, $ref, $alt, $f1, $f2, $f3, $f4;
            $res{"$out.$pop.data"} .= "\n";
        }
    }
}
close VCF;

warn "writing final files\n";
while ( my ($file, $text) = each %res) {
    warn "   $file\n";
    open  O, ">>$file" or die "cannot open $file\n";
    print O $text;
    close O;
}
    
