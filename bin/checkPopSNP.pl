#!/usr/bin/perl

=head1 NAME

checkPopSNP.pl

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

$ARGV[1] or die "usage: PopSNP.var.gz TRIO.var\n";

my $popvar  = shift @ARGV;
my $triovar = shift @ARGV;
my %popvar  = ();
my ($chr, $pos, $rs, $ref, $var, $freq);
my ($f1, $f2, $m1, $m2, $c1, $c2);
my ($fa_p, $ma_p, $ca_p);
my ($fa_f, $ma_f, $ca_f);
my ($fa_r, $ma_r, $ca_r);

warn "loading Population SNPs\n";
open P, "gunzip -c $popvar | " or die "cannot open $popvar\n";
while (<P>) {
    chomp;
    ($chr, $pos, $rs, $ref, $var, $freq) = split (/\t/, $_);
    $chr = "chr$chr";
    $popvar{$chr}{$pos}{'rs'} = $rs;
    my @freq = split (/,/, $freq);
    foreach my $f (@freq) {
        my ($a, $p) = split (/=/, $f);
        $popvar{$chr}{$pos}{$a} = $p;
    }
}
close P;

warn "parsing trio variations\n";
open T, "$triovar" or die "cannot open $triovar\n";
while (<T>) {
    chomp;
    ($chr, $pos, $ref, $f1, $f2, $m1, $m2, $c1, $c2) = split (/\t/, $_);
    $pos++; # 1-based coordinates
    if (defined $popvar{$chr}{$pos}{'rs'}) {
        $fa_f = "$f1/$f2";
        $ma_f = "$m1/$m2";
        $ca_f = "$c1/$c2";
        $fa_r = "$f2/$f1";
        $ma_r = "$m2/$m1";
        $ca_r = "$c2/$c1";
        $fa_p = '-';
        $ma_p = '-';
        $ca_p = '-';
        
        $fa_p = $popvar{$chr}{$pos}{$fa_f} if (defined $popvar{$chr}{$pos}{$fa_f});
        $fa_p = $popvar{$chr}{$pos}{$fa_r} if (defined $popvar{$chr}{$pos}{$fa_r});
        $ma_p = $popvar{$chr}{$pos}{$ma_f} if (defined $popvar{$chr}{$pos}{$ma_f});
        $ma_p = $popvar{$chr}{$pos}{$ma_r} if (defined $popvar{$chr}{$pos}{$ma_r});
        $ca_p = $popvar{$chr}{$pos}{$ca_f} if (defined $popvar{$chr}{$pos}{$ca_f});
        $ca_p = $popvar{$chr}{$pos}{$ca_r} if (defined $popvar{$chr}{$pos}{$ca_r});
        
        $_ .= "$fa_p\t$ma_p\t$ca_p";
    }
    print "$_\n";
}
close T;
