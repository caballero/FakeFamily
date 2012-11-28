#!/usr/bin/perl

=head1 NAME

pedigree.pl

=head1 DESCRIPTION

Simulate a complete pedigree structure using a random selection of founders.

=head1 USAGE

perl pedigree.pl [OPTIONS]

    Parameter          Description                               Default
    -p --pedigree      Pedigree identifier                       [required]
    -s --select        List of pedigrees to use as founders      [random]
    -e --exclude       List of pedigrees to avoid as founders    [none]
    -r --reference     Reference genome version                  [hg19]
    -d --denovo        De novo mutation rate                     [0.0]
    -a --all_ind       Use all possible individuals              [no]
    -o --outdir        Write files in this directory             [localdir]
    -h --help          Print help and exit
    -v --verbose       Verbose mode
    --version          Print version and exit

=head1 EXAMPLES

perl pedigree.pl -p study_id -e unwanted_id

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
use Getopt::Long;
use Pod::Usage;
use List::Util qw/shuffle/;

# Default parameters
my $help     =  undef;         # Print help
my $verbose  =  undef;         # Verbose mode
my $version  =  undef;         # Version call flag
my $pedigree =  undef;
my $select   =  undef;
my $exclude  =  undef;
my $ref      = 'hg19';
my $denovo   =    0.0;
my $demo     =  undef;
my $all_ind  =  undef;
my $outdir   =    '.';

# Main variables
my $our_version  = 0.1;        # Script version number
my $data_table   = '/u5/www/software/gms/current/public/pedigrees/export';
my $reproduction = '/u1/home/jcaballe/FakeFamily/reproduction-cg.pl';
my %exc;
my %sel;
my %ind;
my %fnd;
my %new;
my @bag_fnd_m;
my @bag_fnd_f;
my $nind  = 0;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'version'          => \$version,
    'p|pedigree=s'     => \$pedigree,
    'e|exclude:s'      => \$exclude,
    's|select:s'       => \$select,
    'r|reference:s'    => \$ref,
    'd|denovo:s'       => \$denovo,
    'a|all_ind'        => \$all_ind,
    'o|outdir:s'       => \$outdir,
    'demo'             => \$demo
) or pod2usage(-verbose => 2);

printVersion() if (defined $version);
pod2usage(-verbose => 2) if  (defined $help);
pod2usage(-verbose => 2) if !(defined $pedigree);

if (defined $select) {
    my @sel = split (/,/, $select);
    foreach my $id (@sel) { $sel{$id} = 1; }
}

if (defined $exclude) {
    my @exc = split (/,/, $exclude);
    foreach my $id (@exc) { $exc{$id} = 1; }
}

if ($outdir ne '.') {
    mkdir $outdir unless (-d $outdir);
    die "Problem with output dir $outdir\n" unless (-d $outdir); 
}

loadIndividuals();
swapFounders();
getPathFounders();
linkFounderFiles();
produceOffsprings();

###################################
####   S U B R O U T I N E S   ####
###################################

# printVersion => return version number
sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub loadIndividuals {
    my $file = 'gms_export_individual_pedigrees.txt';
    my %nfnd;
    $nfnd{'female'} = 0;
    $nfnd{  'male'} = 0;
    open F, "$data_table/$file" or die "cannot open $data_table/$file\n";
    while (<F>) {
        chomp;
        next if (m/^Study/); # skip header
        my @line    = split (/\s+/, $_);
        my $ped_id  = $line[1];
        my $ind_id  = $line[3];
        my $sex     = $line[5];
        my $father  = $line[6];
        my $mother  = $line[7];
        my $asm_ref = $line[17];
        next unless (defined $asm_ref);        
        next unless ($asm_ref eq $ref);
        
        if ($ped_id eq $pedigree) {
            if ($father eq 'NULL' or $mother eq 'NULL') {
                $fnd{$ind_id}{'sex'}    = $sex;
                $nfnd{$sex}++;
            }
            else {
                $ind{$ind_id}{'sex'}    = $sex;
                $ind{$ind_id}{'father'} = $father;
                $ind{$ind_id}{'mother'} = $mother;
                $nind++;
            }
            next;
        }
        
        if (defined $exclude) {
            next     if (defined $exc{$ped_id});
        }

        if (defined $select) {
            next unless (defined $sel{$ped_id});
        }
        
        unless (defined $all_ind) {
            next unless ($father eq 'NULL' or $mother eq 'NULL');
        }
        
        if    ($sex eq 'female') {
            push @bag_fnd_f, $ind_id;
        }
        elsif ($sex eq 'male') {
            push @bag_fnd_m, $ind_id;
        }
        
    }
    close F;
    
    if ($nfnd{'female'} > scalar @bag_fnd_f) {
        die "Not enough female founders to swap\n";
    }
    if ($nfnd{  'male'} > scalar @bag_fnd_m) {
        die "Not enough male founders to swap\n";
    }
    
    @bag_fnd_m = shuffle(@bag_fnd_m);
    @bag_fnd_f = shuffle(@bag_fnd_f);
    
    if (defined $verbose) {
        my $fnd  = join ",", keys %fnd;
        my $ind  = join ",", keys %ind;
        my $bagm = join ",", @bag_fnd_m;
        my $bagf = join ",", @bag_fnd_m;
        warn "Pedigree:         $pedigree\n";
        warn "Founders:         $fnd\n";
        warn "Individuals:      $ind\n";
        warn "Elegible males:   $bagm\n";
        warn "Elegible females: $bagf\n";
    }
}

sub swapFounders {
    warn "Swapping founders\n" if (defined $verbose);
    foreach my $ind (keys %fnd) {
        my $sex = $fnd{$ind}{'sex'};
        my $new;
        if ($sex eq 'male') {
            $new = shift @bag_fnd_m;
            $fnd{$ind}{'new'} = $new;
        }
        elsif ($sex eq 'female') {           
            $new = shift @bag_fnd_f;
            $fnd{$ind}{'new'} = $new;
        }
        else {
            die "Error for sex definition for $ind\n";
        }
        
        $new{$new} = $ind;
        warn "$ind <=> $new\n" if (defined $verbose);
    }
}

sub getPathFounders {
    my $file = 'gms_export_assembly_file_pedigrees.txt';
    warn "Finding path to files\n" if (defined $verbose);
    open F, "$data_table/$file" or die "cannot open $data_table/$file\n";
    while (<F>) {
        chomp;
        my $path = undef;
        my $name = undef;
        my @path = ();
        my @line = split (/\s+/, $_);
        my $ind  = $line[3];
        next unless ($line[14] eq 'VAR-ANNOTATION');
        if (defined $new{$ind}) {
            $fnd{ $new{$ind} }{'path'} = $line[15];
        }
        elsif (defined $fnd{$ind}) {
            $path = $line[15];
            @path = split (/\//, $path);
            $name = pop @path;
            $fnd{$ind}{'file'} = $name;
        }
        elsif (defined $ind{$ind}) {
            $path = $line[15];
            @path = split (/\//, $path);
            $name = pop @path;
            $ind{$ind}{'file'} = $name;
        }
    }
    close F;
}

sub linkFounderFiles {
    warn "Linking founders genotypes\n" if (defined $verbose);
    foreach my $ind (keys %fnd) {
        my $file = $fnd{$ind}{'path'};
        my $name = $fnd{$ind}{'file'};
        #system ("cp $file $name");
        system ("ln -s $file $name");
    }
}

sub produceOffsprings {
    my ($ind, $ind_file, $sex, $mother, $mother_file, $father, $father_file);
    foreach $ind (keys %ind) {
        $ind_file = $ind{$ind}{'file'};
        $sex      = $ind{$ind}{'sex'};
        
        next if (-e $ind_file);
        
        $mother = $ind{$ind}{'mother'};
        if    (defined $fnd{$mother}{'file'}) {
            $mother_file = $fnd{$mother}{'file'};
        }
        elsif (defined $ind{$mother}{'file'}) {
            $mother_file = $ind{$mother}{'file'};
        }
        else {
            die "mising mother file for $ind -> $mother\n";
        }
        
        $father = $ind{$ind}{'father'};
        if    (defined $fnd{$father}{'file'}) {
            $father_file = $fnd{$father}{'file'};
        }
        elsif (defined $ind{$father}{'file'}) {
            $father_file = $ind{$father}{'file'};
        }
        else {
            die "mising father file for $ind -> $mother\n";
        }
        
        next unless (-e $mother_file and -e $father_file);
        warn "Generating $ind ($sex), Mother=$mother, Father=$father\n" if (defined $verbose);
        if (defined $demo) {
            system ("touch $ind_file");
        }
        else {
            system ("$reproduction -f $father_file -m $mother_file -s $sex -c $ind_file -d $denovo");
        }
        $nind--;
    }
    
    produceOffsprings() unless ($nind < 1);
}
