#!/usr/bin/perl

=head1 NAME

reproduction.pl

=head1 DESCRIPTION

Do artificial sex and have a baby, you can choose the sex.

=head1 USAGE

perl reproduction.pl [OPTIONS]

OPTIONS:
    Parameter      Description                   Values      Default 
    -f --father    Father data*                  FILE [1]    
    -m --mother    Mother data*                  FILE [1]
    -c --child     Write child data here         FILE [2]    child.out
    -s --sex       Child sex                     M/F  [3]    F
    -a --ave       Average recombination size    SIZE [4]    1Mb
    --min          Minimal recombination size    SIZE [4]    100b
    --max          Maximal recombination size    Size [4]    10Mb
    -g --genome    Reference genome info              [5]    hg19
    -r --reference Reference genome fasta* 
    -n --numsnp    Total number of SNPs          INT         AUTO
    -d --denovo    De novo mutations fraction    0.0-1.0     0.01 (1%)
    -t --hratio    Homo/hetero probility         0.0-1.0     0.5
    
    -h  --help     Print this screen
    -v  --verbose  Verbose mode
    
    * Files can be compressed (gzip/bzip2).
  [1] Input data is a tab-delimited text file with columns (one record 
      per line, sorted by CHROM->POSITION):
      CHROM POSITION REF_ALLELE ALLELE_1 ALLELE_2 ALLELE_INFO
  [2] Output format is tab-delimited text file with columns:
      CHROM POSITION REF_ALLELE ALLELE_1 ALLELE_2 ALLELE_INFO
  [3] Male or Female.
  [4] Sequence sizes can use kb, Mb symbols. Default is b.
  [5] Currently supporting hg19 information (chromosome sizes).

=head1 EXAMPLES

  Simulate a girl:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -r hg19.fa.gz
  
  Simulate a boy, recombination points every 2 Mb
  perl reproduction.pl -f father.var -m mother.var -c boy.var -s M -r hg19.fa.gz --min 2M --max 2M
  
  Simulate a girl, 1M SNPs:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -r hg19.fa.gz -n 1000000
  
  Simulate a girl, no new SNPs:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -d 0
  
  

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

# Global variables
my %size       = ();
my %genome     = ();
my %recom      = ();
my %var        = ();
my @var        = ();
my @allele     = (1, 2);
my %chrom      = ();
my @chrom      = ();
my ($father_h, $mother_h, $genome_h);
my ($chr, $pos, $ref, $al1, $al2, $alt, $info, $allele, $id, $len, $var);
my $nsnp       = 0;
my $newsnp     = 0;
my $hersnp     = 0;
my $num_new_snp= 0;
my @dna        = qw/A C G T/;

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $father     = undef;
my $mother     = undef;
my $child      = 'child.out';
my $sex        = 'F';
my $ave        = 1e6;
my $min        = 100;
my $max        = 1e7;
my $mod        = 'hg19';
my $numsnp     = undef;
my $denovo     = 0.01;
my $genome     = undef;
my $hratio     = 0.5;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'f|father:s'      => \$father,
    'm|mother:s'      => \$mother,
    'c|child:s'       => \$child,
    's|sex:s'         => \$sex,
    'a|ave:s'         => \$ave,
    'min:s'           => \$min,
    'max:s'           => \$max,
    'g|genome:s'      => \$mod,
    'r|reference:s'   => \$genome,
    'n|numsnp:i'      => \$numsnp,
    'd|denovo:s'      => \$denovo,
    't|hratio:s'      => \$hratio
);
pod2usage(-verbose => 2)     if (defined $help);
pod2usage(-verbose => 2) unless (defined $mother and defined $father);

# Calculate number of SNPs to process
$numsnp = calcNumSnp() unless (defined $numsnp);
$newsnp = int($denovo * $numsnp);
$hersnp = $numsnp - $newsnp;
warn "$numsnp will be selected\n" if (defined $verbose);

# Check param
validateParam();

# Load chromosomes sizes
warn "loading chromosome sizes\n" if (defined $verbose);
while (<DATA>) {
    chomp;
    my ($mod, $chr, $len) = split (/:/, $_);
    $size{$mod}{$chr} = $len;
}

# Dirty check of valid models
die "Reference model is't supported $mod\n" unless (defined $size{$mod}{'chr1'});

# Convert size symbols
$ave = convertSymbol($ave);
$min = convertSymbol($min);
$max = convertSymbol($max);

# Meiosis in the Father
warn "doing meiosis in father genome\n" if (defined $verbose);
%recom    = getRecomPoints();
$allele   = $allele[int rand(@allele)]; # initial allele selection, random
$father_h = $father;
$father_h = "gunzip  -c $father | " if ($father =~ m/gz$/);
$father_h = "bunzip2 -c $father | " if ($father =~ m/bz2$/);
open F, "$father_h" or die "cannot open $father\n";
while (<F>) {
    chomp;
    ($chr, $pos, $ref, $al1, $al2, $info) = split (/\t/, $_);
    
    # Sex specific alleles
    next if ($chr eq 'chrX' and $sex eq 'M');
    next if ($chr eq 'chrY' and $sex eq 'F');
    # No recombination in Sex Chromosomes, for now ...
    if ($chr eq 'chrY' or $chr eq 'chrX') {
        $allele = 1;
    }
    else {
        if ($pos > $recom{$chr}->[0]) {
            shift @{ $recom{$chr} };
            # Swap alleles
            if ($allele == 1) { $allele = 2; } else { $allele = 1; }
        }
    }
    $alt = $al1;
    $alt = $al2 if ($allele == 2);    
    $genome{$chr}{$pos}{'father'} = $alt;
    $genome{$chr}{$pos}{'ref'}    = $ref;
    $genome{$chr}{$pos}{'finf'}   = $info unless ($alt eq $ref);
}
close F;

# Meiosis in the Mother
warn "doing meiosis in mother genome\n" if (defined $verbose);
%recom    = getRecomPoints();
$allele   = $allele[int rand(@allele)]; # initial allele selection, random
$mother_h = $mother;
$mother_h = "gunzip  -c $mother | " if ($mother =~ m/gz$/);
$mother_h = "bunzip2 -c $mother | " if ($mother =~ m/bz2$/);
open M, "$mother_h" or die "cannot open $mother\n";
while (<M>) {
    chomp;
    ($chr, $pos, $ref, $al1, $al2, $info) = split (/\t/, $_);
    if ($pos > $recom{$chr}->[0]) {
        shift @{ $recom{$chr} };
        # Swap alleles
        if ($allele == 1) { $allele = 2; } else { $allele = 1; }
    }
    
    $alt = $al1;
    $alt = $al2 if ($allele == 2);    
    $genome{$chr}{$pos}{'mother'} = $alt;
    $genome{$chr}{$pos}{'ref'}    = $ref;
    $genome{$chr}{$pos}{'minf'}   = $info unless ($alt eq $ref);
}
close M;

# Child genome
warn "generating child genome\n" if (defined $verbose);
open STDOUT, ">$child" or die "cannot open $child\n";

foreach $chr (keys %genome) {
    foreach $pos (keys %{ $genome{$chr} }) {
        if ($chr eq 'chrY' and $sex eq 'M') {
            if (defined $genome{$chr}{$pos}{'father'}) {
                $al1 = $genome{$chr}{$pos}{'father'};
                $al2 = '-';
            }
            else {
                next;
            }
        }
        elsif ($chr eq 'chrX' and $sex eq 'M') {
            if (defined $genome{$chr}{$pos}{'mother'}) {
                $al2 = $genome{$chr}{$pos}{'mother'};
                $al1 = '-';
            }
            else {
                next;
            }
        }
        else {
            # Father genotype
            if (defined $genome{$chr}{$pos}{'father'}) {
                $al1 = $genome{$chr}{$pos}{'father'};
            }
            else {
                $al1 = $genome{$chr}{$pos}{'ref'};
            }
        
            # Mother genotype
            if (defined $genome{$chr}{$pos}{'mother'}) {
                $al2 = $genome{$chr}{$pos}{'mother'};
            }
            else {
                $al2 = $genome{$chr}{$pos}{'ref'};
            }
        }
        
        # Reference genotype
        $ref  = $genome{$chr}{$pos}{'ref'};
        
        # Skip reference calls
        next if ($ref eq $al1 and $ref eq $al2);
        
        # Capture RS information
        $info = '';
        if (defined $genome{$chr}{$pos}{'finf'}) {
            if (defined $genome{$chr}{$pos}{'minf'}) {
                $info = $genome{$chr}{$pos}{'finf'} . ',' . $genome{$chr}{$pos}{'minf'};
            }
            else {
                $info = $genome{$chr}{$pos}{'finf'} . ',';
            }
        }
        else {
            if (defined $genome{$chr}{$pos}{'minf'}) {
                $info = ',' .  $genome{$chr}{$pos}{'minf'};
            }
        }  
        
        # Record data
        push @var, "$chr:$pos:$ref:$al1:$al2:$info";
        $nsnp++;
    }
}
undef(%genome); # Free memory
undef(%recom);

# Reduce SNPs if required
#while ($nsnp > $hersnp) {
#    $var = int(rand @var);
#    next unless ($var[$var] =~ m/:/);
#    ($chr,$pos,$ref,$al1,$al2,$info) = split (/:/, $var[$var]);  
#    next if ($al1 eq $al2);   # Prefer to remove heterozygous positions
#    $var[$var] = ''; # deletion
#    $nsnp--;
#}
warn "   $nsnp total parental SNPs\n" if (defined $verbose);

foreach $var (@var) {
    next unless ($var =~ m/:/);
    ($chr, $pos, $ref, $al1, $al2, $info) = split (/:/, $var);
    $var{$chr}{$pos} = "$ref\t$al1\t$al2\t$info";
}
undef(@var);
warn "   $nsnp parental SNPs selected\n" if (defined $verbose);

# Add de novo mutations
if ($nsnp < $numsnp) {
    warn "selecting new SNPs\n" if (defined $verbose);
    
    warn "   loading genome sequences\n" if (defined $verbose);
    %chrom    = ();
    $genome_h = $genome;
    $genome_h = "gunzip  -c $genome | " if ($genome =~ m/gz$/);
    $genome_h = "bunzip2 -c $genome | " if ($genome =~ m/bz2$/);
    open FA, "$genome_h" or die "cannot open $genome!\n";
    $id = undef;
    while (<FA>) {
        chomp;
        if (m/>/) {
            s/>//;
            $id = $_;
        }
        else {
            $chrom{$id} .= $_;
        }
    }
    close FA;

    # chromosome selection depends on size
    foreach $chr (keys %chrom) {
        unless ($chr =~ m/hap|rand|chrM|chrUn/i) {
            for (my $i = 0; $i <= int($size{$mod}{$chr} / 1e6) + 1; $i++) {
                push @chrom, $chr;
            }
        }
    }
    
    warn "   selecting new mutation points\n" if (defined $verbose);
    for (my $i = $nsnp; $i <= $numsnp; $i++) {
        $chr = $chrom[int( rand @chrom)];
        
        if ($chr eq 'chrY' and $sex eq 'F') {
            $i--;
          next;
        }
        
        $len = $size{$mod}{$chr};
	    $pos = int(rand $len);
        
	    if (defined $var{$chr}{$pos}) {
	        $i--;
	      next;
	    }
	
	    $ref = substr ($chrom{$chr}, $pos - 1, 1);
	    if ($ref =~ m/[^ACGT]/) {
	        $i--;
	      next;
	    }

	    my @dna2 = dnaSel($ref, @dna);
	    $al1 = $dna2[int( rand @dna2)];
        $al2 = $al1;

        if ($hratio >= rand) {
             my @dna2 = dnaSel($al1, @dna);
             $al2 = $dna2[int( rand @dna2)];
        }

        if ($chr eq 'chrY' or ($chr eq 'chrX' and $sex eq 'M')) {
            $var{$chr}{$pos} = "$ref\t$al1\t-\t";
        }
        else {
            $var{$chr}{$pos} = "$ref\t$al1\t$al2\t";
        }
        $num_new_snp++;
        $nsnp++;
    }
    undef(%chrom); # free memory
    warn "   $num_new_snp new SNPs added\n" if (defined $verbose);

}

warn "printing SNPs\n" if (defined $verbose);
foreach $chr (sort (keys %var)) {
    foreach $pos (sort {$a<=>$b} (keys %{ $var{$chr} })) {
        print join "\t", $chr, $pos, $var{$chr}{$pos};
        print "\n";
    }
}
warn "done\n" if (defined $verbose);

###########################################################
#            S  U  B  R  O  U  T  I  N  E  S              #
###########################################################
sub convertSymbol {
    my $uni = shift @_;
    my $val = $uni;
    my $fac = 1;
    if ($uni =~ m/(\d+)([Mk])b*/i) {
        $val = $1;
        $fac = uc $2;
        
        if    ($fac eq 'M') { $val *= 1e6; } 
        elsif ($fac eq 'K') { $val *= 1e3; }
        else { die "I cannot recognize $fac in $uni\n"; }
    }
    return $val;
}

sub getRecomPoints {
    my %points = ();
    foreach my $chr (keys %{ $size{$mod} }) {
        my $len = $size{$mod}{$chr};
        my $pos = 0;
        while ($pos <= $len) {
            $pos += getSize();
            push @{ $points{$chr} }, $pos;
        }
        push @{ $points{$chr} }, $len;
    }
    return %points;
}

sub getSize {
    my $size = $ave;
    my $ext  = 0;
    my @op   = ('+', '+', '+', '+', '+', '-',  '-', '-', '-', '-', '++');
    my $op   = $op[int( rand @op)];
    if ($op eq '++') {
        $size += int( rand ($max - $ave));
    }
    elsif ($op eq '+') {
        $size += int( rand ($ave));
    }
    elsif ($op eq '-') {
        $size -= int( rand ($ave));
    }
    else {
        # What am I doing here?
    }
    $size = $min if ($size < $min);
    $size = $max if ($size > $max);
    return $size;
}

sub calcNumSnp {
    my $num   = 3e6;
    my $noise = int( rand( $num * 0.1));
    #my @ops   = ( '+', '-');
    my @ops   = ( '+', '+');
    my $op    = $ops[int( rand @ops)];
    if    ($op eq '+') { $num += $noise; }
    elsif ($op eq '-') { $num -= $noise; }
    return $num;
}

sub dnaSel {
    my $rem = shift @_;
    my @res = ();
    foreach my $n (@_) {
        push @res, $n unless ($n eq $rem);
    }
    return @res;
}

sub validateParam {
    die "father file not found: $father\nKilled\n"    unless (-e $father);
    die "mother file not found: $mother\nKilled\n"    unless (-e $mother);
    die "numsnp isn't numeric:   $numsnp\nKilled\n"   unless ($numsnp =~ m/\d+/ and $numsnp > 1);
    die "denovo must be 0.0-1.0: $denovo\nKilled\n"   unless ($denovo =~ m/\d+/ and $denovo >= 0 and $denovo <= 1);
    die "hratio must be 0.0-1.0: $hratio\nKilled\n"   unless ($hratio =~ m/\d+/ and $hratio >= 0 and $hratio <= 1);
    die "sex must be F or M:     $sex\nKilled\n"      unless ($sex eq 'F' or $sex eq 'M');
}

__END__
hg19:chr7:159138663
hg19:chr20:63025520
hg19:chr22:51304566
hg19:chr14:107349540
hg19:chrY:59373566
hg19:chr8:146364022
hg19:chr19:59128983
hg19:chr1:249250621
hg19:chr6:171115067
hg19:chr11:135006516
hg19:chr17:81195210
hg19:chr21:48129895
hg19:chr16:90354753
hg19:chr18:78077248
hg19:chr3:198022430
hg19:chr12:133851895
hg19:chr15:102531392
hg19:chrX:155270560
hg19:chr4:191154276
hg19:chrM:16571
hg19:chr2:243199373
hg19:chr9:141213431
hg19:chr13:115169878
hg19:chr10:135534747
hg19:chr5:180915260
