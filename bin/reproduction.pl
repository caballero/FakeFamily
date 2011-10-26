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
    -p --points    Total recombination points    INT         
    
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
  [5] Currently supporting hg19/hg18 information (chromosome sizes).

=head1 EXAMPLES

  Simulate a girl:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -r hg19.fa.gz
  
  Simulate a boy, recombination points every 2 Mb
  perl reproduction.pl -f father.var -m mother.var -c boy.var -s M -r hg19.fa.gz --min 2M --max 2M
  
  Simulate a girl, 1M SNPs:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -r hg19.fa.gz -n 1000000
  
  Simulate a girl, no new SNPs:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -d 0
  
  Simulate a girl, no new SNPs, max 20 recombination points:
  perl reproduction.pl -f father.var -m mother.var -c girl.var -d 0 -p 20

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
use List::Util qw/shuffle/;
# sorry I need to hardcode the HapMap lib, in a future we'll share the lib and data
if (-e '/proj/famgen/lib/gglusman/HapMap.pm') { #bobama
    use lib "/proj/famgen/lib/gglusman";
}
if (-e '/proj/famgen/lib/gglusman/HapMap.pm') { #osiris
    use lib "/net/gestalt/system/utils";
}
use HapMap;

# Global variables
my %size       = ();
my %centro     = ();
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
my $points     = undef;

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
    't|hratio:s'      => \$hratio,
    'p|points:s'      => \$points
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
    next if (/#/);
    my ($mod, $chr, $len, $cen) = split (/:/, $_);
    $size{$mod}{$chr} = $len;
    $centro{$mod}{$chr} = $cen;
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
    next if ($chr eq 'chrM');
    
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
    unless ($chr =~ m/chrM|chrY|chrX/) {
        if ($pos > $recom{$chr}->[0]) {
            shift @{ $recom{$chr} };
            # Swap alleles
            if ($allele == 1) { $allele = 2; } else { $allele = 1; }
        }
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
        elsif ($chr eq 'chrM') {
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
while ($nsnp > $hersnp) {
    $var = int(rand @var);
    next unless ($var[$var] =~ m/:/);
#    ($chr,$pos,$ref,$al1,$al2,$info) = split (/:/, $var[$var]);  
#    next if ($al1 eq $al2);   # Prefer to remove heterozygous positions
    $var[$var] = ''; # deletion
    $nsnp--;
}
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
    my $num    = 0;
    foreach my $chr (keys %{ $size{$mod} }) {
        next if ($chr =~ m/chrM|chrY|chrX/);
        my $len = $size{$mod}{$chr};
        my @pos = ();
        my $pos = 0;
        while ($pos <= $len) {
            $pos += getSize();
            push @pos, $pos;
        }
        my @pos_fil = filterPointsHapMap($chr, @pos);
        @{ $points{$chr} } = @pos_fil;
        $num += $#pos_fil + 1;
    }
    
    warn "recom $num, wanted $points\n" if (defined $verbose);
    if (defined $points) {
        if ($num > $points) {
            # reduce the recombination points, random deletions
            my @pos = ();
            foreach my $chr (keys %points) {
                foreach my $pos (@{ $points{$chr} }) {
                    push @pos, "$chr:$pos";
                }
            }
            my $pp = join (",", @pos);
            warn "original points: $pp\n" if (defined $verbose);
            @pos    = shuffle(@pos);
            $pp = join (",", @pos);
            warn "shuffle points: $pp\n" if (defined $verbose);
            %points = ();
            for (my $i = 0; $i < $points; $i++) {
                my ($chr, $pos) = split (/:/, $pos[$i]);
                push @{ $points{$chr} }, $pos;
            }
        }
    }
    # add tails to the data
    foreach my $chr (keys %{ $size{$mod} }) {
        next if ($chr =~ m/chrM|chrY|chrX/);
        my $len = $size{$mod}{$chr};
        push @{ $points{$chr} }, $len;
        my $pts = length @{ $points{$chr} };
        warn "    $pts recombination points in $chr\n" if (defined $verbose);
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

sub filterPoints {
    my $chr    = shift @_;
    my $centro = $centro{$mod}{$chr};
    my $size   = $size{$mod}{$chr};
    my $left   = $centro;
    my $right  = $size - $centro;
    my @pos = ();
	warn "    original recombination points $#_\n" if (defined $verbose);
    foreach my $pos (@_) {
        my $dice = rand();
        my $prob = 0;
        if ($pos < $centro) {
            $prob = 1 - ($pos / $left);
        }
        elsif ($pos > $centro) {
            $prob = $pos / $right;
        }
        
        push @pos, $pos if ($dice < $prob);
    }
	warn "    filtered recombination points $#pos\n" if (defined $verbose);
    return @pos;
}

sub filterPointsHapMap {
    my $chr  = shift @_;
	my $hm   = new HapMap($mod);
    my @pos  = ();
    my $lim  = 7; # minimal genetic distance to consider a recombination
	warn "    $chr original recombination points $#_\n" if (defined $verbose);
	my $last = 1;
    foreach my $pos (@_) {
        my $gendis = $hm->genetic_distance($chr, $last, $pos);
		if ($gendis > $lim) { 
			push @pos, $pos;
			$last    = $pos;
		}
    }
	warn "    $chr filtered recombination points $#pos\n" if (defined $verbose);
    return @pos;
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
#build:chr:size:centromere
hg18:chr1:247249719:124300000
hg18:chr2:242951149:93300000
hg18:chr3:199501827:91700000
hg18:chr4:191273063:50700000
hg18:chr5:180857866:47700000
hg18:chr6:170899992:60500000
hg18:chr7:158821424:59100000
hg18:chr8:146274826:45200000
hg18:chr9:140273252:51800000
hg18:chr10:135374737:40300000
hg18:chr11:134452384:52900000
hg18:chr12:132349534:35400000
hg18:chr13:114142980:16000000
hg18:chr14:106368585:15600000
hg18:chr15:100338915:17000000
hg18:chr16:88827254:38200000
hg18:chr17:78774742:22200000
hg18:chr18:76117153:16100000
hg18:chr19:63811651:28500000
hg18:chr20:62435964:27100000
hg18:chr21:46944323:12300000
hg18:chr22:49691432:11800000
hg18:chrX:154913754:59500000
hg18:chrY:57772954:11300000
hg18:chrM:16571:16571
hg19:chr1:249250621:125000000
hg19:chr2:243199373:93300000
hg19:chr3:198022430:91000000
hg19:chr4:191154276:50400000
hg19:chr5:180915260:48400000
hg19:chr6:171115067:61000000
hg19:chr7:159138663:59900000
hg19:chr8:146364022:45600000
hg19:chr9:141213431:49000000
hg19:chr10:135534747:40200000
hg19:chr11:135006516:53700000
hg19:chr12:133851895:35800000
hg19:chr13:115169878:17900000
hg19:chr14:107349540:17600000
hg19:chr15:102531392:19000000
hg19:chr16:90354753:36600000
hg19:chr17:81195210:24000000
hg19:chr18:78077248:17200000
hg19:chr19:59128983:26500000
hg19:chr20:63025520:27500000
hg19:chr21:48129895:13200000
hg19:chr22:51304566:14700000
hg19:chrX:155270560:60600000
hg19:chrY:59373566:12500000
hg19:chrM:16571:16571
