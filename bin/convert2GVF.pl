#!/usr/bin/perl

=head1 NAME

convert2GVF.pl 

=head1 DESCRIPTION

Parse a personal genome file (SNPs) and produce the equivalent GVF file.

=head1 USAGE

OPTIONS
    Parameter        Description                Value      Default
    -i --input       Input                      File       STDIN
    -o --output      Output                     File       STDOUT
    -r --reference   Genome reference version   Name       hg19
    -q --qual        Quality score              INT        .
    -l --label       Label for source           Name       fake
    -h --help        Print this screen
    -v --verbose     Verbose mode

=head1 EXAMPLES

    

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
my $input      = undef;
my $output     = undef;
my $ref        = 'hg19';
my $qual       = '.';
my $label      = 'fake';
# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'i|input:s'       => \$input,
    'o|output:s'      => \$output,
    'r|reference:s'   => \$ref, 
    'q|qual:i'        => \$qual,
    'l|label:s'       => \$label
);
pod2usage(-verbose => 2) if (defined $help);

if (defined $input) {
    open STDIN, "$input" or die "Cannot read file $input\n";
}

if (defined $output) {
    open STDOUT, ">$output" or die "Cannot write file $output\n";
}

my $num = 0;
my $header = getHeader($ref);
print "$header\n";
while (<>) {
    chomp;
    $num++;
    my ($chr, $pos, $r1, $a1, $a2, $info) = split (/\t/, $_);
    my $gtype = 'heterozygous';
    my $gtype = 'homozygous' if ($a1 eq $a2);
    my $lab   = $label;
    if ($info =~ /rs\d+/) {
        $lab = $info;
        $lab =~ s/,$//;
        $lab =~ s/^,//;
    }
    
    my $tag = "ID=SNV$num;Variant_seq=$a1,$a2;Reference_seq=$r1;Genotype=$gtype;";
    print join "\t", $chr, $lab, 'SNV', $pos, $pos, $qual, '+', '.', $tag;
    print "\n";
}

sub getHeader {
    my $ref = shift @_;
    my $head = "##gff-version 3\n##gvf-version 1.02\n";
    
    if ($read eq 'hg19') {
        $head .= <<_last
##sequence-region chr1  1 247249719
##sequence-region chr2  1 242951149
##sequence-region chr3  1 199501827
##sequence-region chr4  1 191273063
##sequence-region chr5  1 180857866
##sequence-region chr6  1 170899992
##sequence-region chr7  1 158821424
##sequence-region chr8  1 146274826
##sequence-region chr9  1 140273252
##sequence-region chr10 1 135374737
##sequence-region chr11 1 134452384
##sequence-region chr12 1 132349534
##sequence-region chr13 1 114142980
##sequence-region chr14 1 106368585
##sequence-region chr15 1 100338915
##sequence-region chr16 1 88827254
##sequence-region chr17 1 78774742
##sequence-region chr18 1 76117153
##sequence-region chr19 1 63811651
##sequence-region chr20 1 62435964
##sequence-region chr21 1 46944323
##sequence-region chr22 1 49691432
##sequence-region chrX  1 154913754
##sequence-region chrY  1 57772954
##sequence-region chrM  1 16571
_last
;   }
    elsif ($read eq 'GRCh37') {
        $head .= <<_last
##sequence-region chr1  1 247249719
##sequence-region chr2  1 242951149
##sequence-region chr3  1 199501827
##sequence-region chr4  1 191273063
##sequence-region chr5  1 180857866
##sequence-region chr6  1 170899992
##sequence-region chr7  1 158821424
##sequence-region chr8  1 146274826
##sequence-region chr9  1 140273252
##sequence-region chr10 1 135374737
##sequence-region chr11 1 134452384
##sequence-region chr12 1 132349534
##sequence-region chr13 1 114142980
##sequence-region chr14 1 106368585
##sequence-region chr15 1 100338915
##sequence-region chr16 1 88827254
##sequence-region chr17 1 78774742
##sequence-region chr18 1 76117153
##sequence-region chr19 1 63811651
##sequence-region chr20 1 62435964
##sequence-region chr21 1 46944323
##sequence-region chr22 1 49691432
##sequence-region chrX  1 154913754
##sequence-region chrY  1 57772954
##sequence-region chrM  1 16571
_last
;   }
    elsif ($read eq 'hg18') {
        $head .= <<_last
##sequence-region chr1  1 247249719
##sequence-region chr2  1 242951149
##sequence-region chr3  1 199501827
##sequence-region chr4  1 191273063
##sequence-region chr5  1 180857866
##sequence-region chr6  1 170899992
##sequence-region chr7  1 158821424
##sequence-region chr8  1 146274826
##sequence-region chr9  1 140273252
##sequence-region chr10 1 135374737
##sequence-region chr11 1 134452384
##sequence-region chr12 1 132349534
##sequence-region chr13 1 114142980
##sequence-region chr14 1 106368585
##sequence-region chr15 1 100338915
##sequence-region chr16 1 88827254
##sequence-region chr17 1 78774742
##sequence-region chr18 1 76117153
##sequence-region chr19 1 63811651
##sequence-region chr20 1 62435964
##sequence-region chr21 1 46944323
##sequence-region chr22 1 49691432
##sequence-region chrX  1 154913754
##sequence-region chrY  1 57772954
##sequence-region chrM  1 16571
_last
;   }
    elsif ($read eq 'NCBI B36') {
        $head .= <<_last
##sequence-region chr1  1 247249719
##sequence-region chr2  1 242951149
##sequence-region chr3  1 199501827
##sequence-region chr4  1 191273063
##sequence-region chr5  1 180857866
##sequence-region chr6  1 170899992
##sequence-region chr7  1 158821424
##sequence-region chr8  1 146274826
##sequence-region chr9  1 140273252
##sequence-region chr10 1 135374737
##sequence-region chr11 1 134452384
##sequence-region chr12 1 132349534
##sequence-region chr13 1 114142980
##sequence-region chr14 1 106368585
##sequence-region chr15 1 100338915
##sequence-region chr16 1 88827254
##sequence-region chr17 1 78774742
##sequence-region chr18 1 76117153
##sequence-region chr19 1 63811651
##sequence-region chr20 1 62435964
##sequence-region chr21 1 46944323
##sequence-region chr22 1 49691432
##sequence-region chrX  1 154913754
##sequence-region chrY  1 57772954
##sequence-region chrM  1 16571
_last
;   }
    else {
        die "I don't have genome information for $ref, aborting\n";   
    }
    
    return $head;
}
