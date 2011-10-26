#!/usr/bin/perl

=head1 NAME

noise.pl

=head1 DESCRIPTION

Add noise (errors in sequencing) in a variation file.

=head1 USAGE

perl noise.pl [OPTIONS]

OPTIONS
   Parameter                Description              Value       Default
    -i --input              Input  VAR file*         FILE        
    -o --output             Output VAR file          FILE
    -g --genome             Genome fasta file*       FILE
    -s --sex                Sex to use               M/F***

   Error rates:
    --ur --undercall_ref    Undercall reference      Float       1e-3
    --un --undercall_noref  Undercall no reference   Float       2e-4
    --or --overcall_ref     Overcall reference       Float       2e-6
    --on --overcall_noref   Overcall no reference    Float       2e-4
    
    --nonew       No new mutations**
    
    -h --help     Print this screen
    -v --verbose  Verbose mode on
    
     * Files can be compressed gzip|bzip2
    ** If used, a reference genome isn't required
   *** Male or Female 
    
=head1 EXAMPLES

    # Typical use
    perl noise.pl -g hg19.fa.gz -i genome.var -o mod_genome.var
    
    # Increase undercall no reference error rate to 1e-2
    perl noise.pl -i genome.var -o mod_genome.var --un 0.01 -g hg19.fa.gz 

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
my $genome     = undef;
my $nonew      = undef;
my $sex        = undef;
my $uc_ref     = 1e-3;
my $uc_noref   = 2e-4;
my $oc_ref     = 2e-6;
my $oc_noref   = 2e-4;

# Fetch options
GetOptions(
    'h|help'             => \$help,
    'v|verbose'          => \$verbose,
    'i|input:s'          => \$input,
    'o|output:s'         => \$output,
    'g|genome:s'         => \$genome,
    'n|nonew'            => \$nonew,
    's|sex:s'            => \$sex,
    'uc|undercall_ref'   => \$uc_ref,
    'un|undercall_noref' => \$uc_noref,
    'oc|overcall_ref'    => \$oc_ref,
    'on|overcall_noref'  => \$oc_noref,
);
pod2usage(-verbose => 2) if (defined $help);
pod2usage(-verbose => 2) unless (defined $input and defined $output and defined $sex);

if (defined $output) {
    open STDOUT, ">$output" or die "cannot open $output\n";
}

# Global variables
my ($file_h);
my ($chr, $pos, $ref, $al1, $al2, $info, $err, $var);
my @var        = ();
my @chrom      = ();
my @dna        = qw/A C G T/;
my @dna2       = ();
my $numerr     = 0;
my @error      = ();
my %genome     = ();
my %size       = ();
my %var        = ();
my %change     = ();

# Reading variation points 
warn "reading variation data\n" if (defined $verbose);
$file_h = $input;
$file_h = "gunzip  -c $input | " if ($input =~ m/gz$/);
$file_h = "bunzip2 -c $input | " if ($input =~ m/bz2$/);
open IN, "$file_h" or die "cannot open $input\n";
@var = <IN>;
close IN;

# loading genome asequences
loadGenome() unless (defined $nonew);

# Computing number or errors
$numerr = calcErrors($#var);
warn "perfoming $numerr changes\n" if (defined $verbose);

# Doing changes
for (my $i = 0; $i <= $numerr; $i++) {
    $err = $error[int(rand @error)];
    $var = int(rand @var);
    if (defined $change{$var}) {
        $i--;
        next;
    }
    $change{$var} = 1;
    
    if ($err eq 'del') {
        $var[$var] = undef;
    }
    elsif ($err eq 'sus') {
        ($chr, $pos, $ref, $al1, $al2, $info) = split (/\t/, $var[$var]);
        if ($sex eq 'F' or $chr =~ m/chr\d+/) {
            # random selection of one allele
            if (0.5 > rand) {
                @dna2 = dnaSel($al1, @dna);
                $al1  = $dna[int(rand @dna2)];
            } 
            else {
                @dna2 = dnaSel($al2, @dna);
                $al2  = $dna[int(rand @dna2)];
            }
        
            if ($al1 eq $ref and $al2 eq $ref) { # Homozygous to reference
                $var[$var] = undef;
            }
            else {
                $var[$var] = join ("\t", $chr, $pos, $ref, $al1, $al2, $info);
            }
        }
        else { # it's a Male genome and it's chrX or chrY
            if ($al1 =~ m/[ACGT]/) {
                @dna2 = dnaSel($al1, @dna);
                $al1  = $dna[int(rand @dna2)];
            }
            else {
                @dna2 = dnaSel($al2, @dna);
                $al2  = $dna[int(rand @dna2)];
            }
                
            if ($al1 eq $ref or $al2 eq $ref) { # Homozygous to reference
                $var[$var] = undef;
            }
            else {
                $var[$var] = join ("\t", $chr, $pos, $ref, $al1, $al2, $info);
            }
        }   
    }
    elsif ($err eq 'new') {
        my $ok = 0;
        while ($ok == 0) {
            $chr  = $chrom[int(rand @chrom)];
            next if ($chr eq 'chrY' and $sex eq 'F');
            $pos  = int(rand $size{$chr});
            $ref  = substr($genome{$chr}, $pos, 1);
            next unless ($ref =~ m/[ACGT]/);
            @dna2 = dnaSel($ref, @dna);
            $al1  = $dna2[int(rand @dna2)];
            if (0.5 > rand) {
                $al2 = $al1;
            }
            else {
                @dna2 = dnaSel($al1, @dna);
                $al2  = $dna2[int(rand @dna2)];
            }
            if ($chr =~ m/chr[XY]/ and $sex eq 'M') {
                push @var, "$chr\t$pos\t$ref\t$al1\t-\n";
            }
            else {
                push @var, "$chr\t$pos\t$ref\t$al1\t$al2\n";
            }
            $ok = 1;
        }
    }
}

warn "sorting new variations\n" if (defined $verbose);
foreach $var (@var) {
    next unless (defined $var);
    ($chr, $pos) = split (/\t/, $var);
    $var{$chr}{$pos} = $var;
}
foreach $chr (sort(keys %var)) {
    foreach $pos (sort {$a<=>$b} keys %{ $var{$chr} }) {
        print $var{$chr}{$pos};
    }
}
warn "done\n" if (defined $verbose);

###########################################################
#            S  U  B  R  O  U  T  I  N  E  S              #
###########################################################
sub dnaSel {
    my $rem = shift @_;
    my @res = ();
    foreach my $n (@_) {
        push @res, $n unless ($n eq $rem);
    }
    return @res;
}

sub calcErrors {
    my $tot = shift @_;
    my $ur  = int($uc_ref   * $tot);
    my $un  = int($uc_noref * $tot);
    my $or  = int($oc_ref   * $tot);
    my $on  = int($oc_noref * $tot);
    my $sum = $ur + $or + $un + $on;
    return $sum;
}

sub loadGenome {
    warn "loading genome sequences\n" if (defined $verbose);
    $file_h = $genome;
    $file_h = "gunzip  -c $genome | " if ($genome =~ m/gz$/);
    $file_h = "bunzip2 -c $genome | " if ($genome =~ m/bz2$/);
    open F, "$file_h" or die "cannot open $genome\n";
    while (<F>) {
        chomp;
        if (m/>(.+)/) {
            $chr = $1;
            push @chrom, $chr;
        }
        else {
            $genome{$chr} .= $_;
            $size{$chr} += length $_;
        }
    }
    close F;
}
