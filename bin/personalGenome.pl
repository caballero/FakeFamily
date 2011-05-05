#!/usr/bin/perl
$| = 1;

=head1 NAME

personalGenome.pl

=head1 DESCRIPTION

Simulate an individual genome file. Our strategy is to use common SNPs 
reported,selecting by allele frequencies and add some de novo mutations.
The final product is a simple table with all the mutations added.

=head1 USAGE

perl createPersonalGenome.pl [OPTIONS]

OPTIONS
    Parameter       Description                 Values       Default 
    -g  --genome    Genome file*                FILE   [1]   
    -p  --popvar    Population variation data*  FILE   [2]
    -d  --dbsnp     dbSNP variation data*       FILE   [3]
    -n  --numsnp    Total SNPS to add           INT          AUTO
    -r  --hratio    Homo/hereto probability     0.0-1.0      0.5
    -s  --sex       Sex of the individual       M/F    [4]   F
    -o  --output    Write output table here     FILE   [5]   person.out
    -w  --newsnp    Total de-novo SNPS to add   0.0-1.0      0.01 (1%)
    
    -k  --keep_ref  Report reference alleles     
        --no_ref    Don't use a reference genome       [6]
        --no_popvar Don't use population  data         [7]
        --no_dbsnp  Don't use dbSNP variation data     [7]
    
    -h  --help         Print this screen
    -v  --verbose      Verbose mode
    
    * Files can be compressed (gzip/bzip2).
  [1] Fasta format.
  [2] Variation data is a tab-delimited text file with columns (one record 
      per line):
      CHROM POSITION REF_ALLELE ALT_ALLELE ALLELE_FREQ
  [3] Variation data is a tab-delimited text file with columns (one record 
      per line):
      CHROM POSITION REF_ALLELE ALT_ALLELE ALLELE_FREQ
  [4] Male or Female.
  [5] Output format is tab-delimited text file with columns:
      CHROM POSITION REF_ALLELE ALLELE_1 ALLELE_2 ALLELE_INFO
  [6] If you don't use a reference, de novo SNPs aren't created.
  [7] If you don't use dbSNP and population variation data, all SNPs will 
      be random.

=head1 EXAMPLES

  # Create a JPT genome, Female
  perl personalGenome.pl -g hg19.fa.gz -d dbSNP132.var.gz -p JPT.var.gz -s F

  # Create a CEU genome, Male, no de-novo variants
  perl personalGenome.pl -p CEU.var.gz -d dbSNP132.var.gz -s M --no_ref
  
  # Create a genome, Male, all random SNPs homozygous
  perl personalGenome.pl -g hg19.fa.gz -s M -r 0 --no_pop --no_dbsnp
  
  # Create a CHB genome, Female, only 1.5 million SNPs
  perl personalGenome.pl -g hg19.fa.gz -d dbSNP132.var.gz -p CHB.var.gz -n 1500000
  
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
my @dna        = qw/A C G T/;
my @chrom      = ();
my %chrom      = ();
my @var        = ();
my %var        = ();
my %known      = ();
my %len        = ();
my $nsnp       =  0;
my $num_pop_snp=  0;
my $num_db_snp =  0;
my $num_new_snp=  0;

my ($chr, $pos, $rs, $alt, $al1, $al2, $ref, $af, $gmaf, $id, $len, $var);
my ($genome_h, $popvar_h, $dbsnp_h);

# Parameters initialization
my $help       = undef;
my $verbose    = undef;
my $genome     = undef;
my $popvar     = undef;
my $dbsnp      = undef;
my $numsnp     = undef;
my $numsnpk    = undef;
my $newsnp     = 0.01;
my $hratio     = 0.5;
my $sex        = 'F';
my $output     = 'person.out';
my $noref      = undef;
my $nopop      = undef;
my $nodbsnp    = undef;
my $keepref    = undef;

# Fetch options
GetOptions(
    'h|help'          => \$help,
    'v|verbose'       => \$verbose,
    'g|genome:s'      => \$genome,
    'p|popvar:s'      => \$popvar,
    'd|dbsnp:s'       => \$dbsnp,
    'n|numsnp:i'      => \$numsnp,
    'w|newsnp:s'      => \$newsnp,
    'r|hratio:s'      => \$hratio,
    's|sex:s'         => \$sex,
    'o|output:s'      => \$output,
    'k|keepref'       => \$keepref,
    'no_ref'          => \$noref,
    'no_popvar'       => \$nopop,
    'no_dbsnp'        => \$nodbsnp
);
pod2usage(-verbose => 2)     if (defined $help);
pod2usage(-verbose => 2) unless (defined $popvar or defined $genome or defined $dbsnp);
$numsnp = calcNumSnp() unless (defined $numsnp);
validateParam();

# Output file
if (defined $output) {
    open STDOUT, ">$output" or die "cannot write output file $output\n";
}

# Computing new and known SNPs expected
warn "$numsnp SNPs to be added\n" if (defined $verbose);
if (defined $noref) {
    $numsnpk = $numsnp;
}
else {
    $numsnpk = $numsnp - int($newsnp * $numsnp);
}

# Population variation SNPs selection
unless (defined $nopop) {
    warn "selecting population SNPs\n" if (defined $verbose);
    warn "   loading known SNPs records\n" if (defined $verbose);
    $popvar_h = $popvar;
    $popvar_h = "gunzip  -c $popvar | " if ($popvar =~ m/gz$/);
    $popvar_h = "bunzip2 -c $popvar | " if ($popvar =~ m/bz2$/);
    open PV, "$popvar_h" or die "cannot open $popvar!\n";
    @var = <PV>;
    close PV;    
    shuffle(\@var);
    
    warn "   selecting SNPs\n" if (defined $verbose);
    while ($var = shift @var) {
        next if ($var =~ m/^#/);
        chomp $var;
        ($chr, $pos, $rs, $ref, $alt, $af) = split (/\t/, $var);
        $chr     = "chr$chr";
        next if (length $ref > 1); # no indels for now
        next if ($chr eq 'chrMT');
        next if ($chr eq 'chrY' and $sex eq 'F');
        my @alt  = split (/,/, $alt);
        my @af   = split (/,/, $af);
        my $geno = "$ref/$ref";
        my $facc = 0;
        my $dice = rand;
        foreach my $alle (@af) {
            my ($g, $f) = split (/=/, $alle);
            $facc += $f;
            if ($dice < $facc) { 
                $geno = $g;
            }
            else {
                last;
            }
        }
        
        unless (defined $keepref) {
            next if ($geno eq "$ref/$ref");
        }
        
        ($al1, $al2) = split (/\//, $geno);
        next if (length $al1 > 1 or length $al2 > 1); # only SNPs
        
        if ($chr eq 'chrY' or ($chr eq 'chrX' and $sex eq 'M')) {
            my $ale = $al1;
            $ale = $al2 if ($al1 eq $ref);
            $var{$chr}{$pos} = "$ref\t$ale\t-\t$rs";
        }
        else {
            $var{$chr}{$pos} = "$ref\t$al1\t$al2\t$rs";
        }
        $nsnp++;
        $num_pop_snp++;
        last if ($nsnp >= $numsnpk);
    }
    warn "   $num_pop_snp population SNPs selected\n" if (defined $verbose);
    undef(@var); # free memory
}

# dbSNP variation SNPs selection
unless (defined $nodbsnp or $nsnp >= $numsnpk) {
    warn "selecting SNPs from dbSNPs\n" if (defined $verbose);
    warn "   loading dbSNPs records\n" if (defined $verbose);
    $dbsnp_h = $dbsnp;
    $dbsnp_h = "gunzip  -c $dbsnp | " if ($dbsnp =~ m/gz$/);
    $dbsnp_h = "bunzip2 -c $dbsnp | " if ($dbsnp =~ m/bz2$/);
    open DB, "$dbsnp_h" or die "cannot open $dbsnp!\n";
    # Take random 1/3 of the dbSNP to save memory
    while (<DB>) {
        next if (m/^[MT|PAR]/); # exclude MT and PAR chromosomes
        push @var, $_ if (0.3333 > rand);   
    }
    close DB;
    warn "   dbSNP subsample: $#var\n" if (defined $verbose);
    shuffle(\@var);
    
    warn "   selecting SNPs\n" if (defined $verbose);
    while ($var = shift @var) {
        next if ($var =~ m/^#/);
        chomp $var;
        ($chr, $pos, $rs, $ref, $alt, $af) = split (/\t/, $var);
        $chr     = "chr$chr";
        next if (defined $var{$chr}{$pos});
        next if (length $ref > 1); # no indels for now
        next if ($chr eq 'chrMT');
        next if ($chr eq 'chrY' and $sex eq 'F');
        $af = rand if ($af eq '-'); # Umm no value?, not problem, random
        next unless ($af >= rand);
        my @alt  = split (/,/, $alt);
        
        $al1     = $alt[int (rand @alt)];
        $al2     = $al1;
        
        if ($hratio > rand) {
            my @dna2 = dnaSel($al1, @dna);
            $al2 = $dna2[int( rand @dna2)];
        }
        
        next if (length $al1 > 1 or length $al2 > 1); # only SNPs
        
        if ($chr eq 'chrY' or ($chr eq 'chrX' and $sex eq 'M')) {
            my $ale = $al1;
            $ale = $al2 if ($al1 eq $ref);
            $var{$chr}{$pos} = "$ref\t$ale\t-\t$rs";
        }
        else {
            $var{$chr}{$pos} = "$ref\t$al1\t$al2\t$rs";
        }
        $nsnp++;
        $num_db_snp++;
        last if ($nsnp >= $numsnpk);
    }
    warn "   $num_db_snp known SNPs selected\n" if (defined $verbose);
    undef(@var); # free memory
}

# de novo SNPs selection
unless (defined $noref or $nsnp >= $numsnp) {
    warn "selecting new SNPs\n" if (defined $verbose);
    
    warn "   loading genome sequences\n" if (defined $verbose);
    %chrom    = ();
    $genome_h = $genome;
    $genome_h = "gunzip  -c $genome | " if ($genome =~ m/gz$/);
    $genome_h = "bunzip2 -c $genome | " if ($genome =~ m/bz2$/);
    open FA, "$genome_h" or die "cannot open $genome : $genome_h!\n";
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

    foreach $chr (keys %chrom) {
        unless ($chr =~ m/hap|rand|chrM|chrUn/i) {
            $len{$chr} = length $chrom{$chr};
            for (my $i = 0; $i <= int($len{$chr} / 1e6) + 1; $i++) {
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
        
        $len = $len{$chr};
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
 warn "$nsnp total SNPs\n" if (defined $verbose);

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
sub calcNumSnp {
    my $num   = 3e6;
    my $noise = int( rand( $num * 0.1));
    my @ops   = ( '+', '-');
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

# Fisher-Yates shuffle algorithm
sub shuffle {
    my $deck = shift; # $deck is a reference to an array
    return unless @$deck; # must not be empty!
    my $i = @$deck;
    while (--$i) {
        my $j = int rand ($i+1);
        @$deck[$i,$j] = @$deck[$j,$i];
    }
}

sub validateParam {
  if (defined $genome) {
    die "Genome not found: $genome\nKilled\n"         unless (-e $genome);
  }
  if (defined $popvar) {
    die "Pop. variation not found: $popvar\nKilled\n" unless (-e $popvar);
  }
  if (defined $dbsnp) {
    die "dbDNP data not found: $dbsnp\nKilled\n"      unless (-e $dbsnp);
  }
    die "numsnp isn't numeric:   $numsnp\nKilled\n"   unless ($numsnp =~ m/\d+/ and $numsnp > 1);
    die "newsnp must be 0.0-1.0: $newsnp\nKilled\n"   unless ($newsnp =~ m/\d\.\d+/ and $newsnp >= 0 and $newsnp <= 1);
    die "hratio must be 0.0-1.0: $hratio\nKilled\n"   unless ($hratio =~ m/\d\.\d+/ and $hratio >= 0 and $hratio <= 1);
    die "sex must be F or M:     $sex\nKilled\n"      unless ($sex eq 'F' or $sex eq 'M');
}
