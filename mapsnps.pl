#!/usr/bin/env perl
use warnings;
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=pod

=head1 NAME

mapsnps.pl - Map genomic SNPs to known human genes / proteins

=head1 SYNOPSIS

mapsnsp.pl [options] infile [>outfile]

where options are:

  -c          annotate coding SNPs only

  -m          annotate missense coding SNPs only

  -n          skip alleles without a matching reference nt

  -y file     save PolyPhen-2 queries for missense SNPs to file

  -g hgXX     genome assembly version used for annotations; default
              is hg19 (GRCh37), also supported is hg18 (NCBI36)

  -i hgYY     input SNP chromosome coordinates assembly version; default
              is to assume the same assembly as used for the annotations,
              otherwise input coordinates will be transformed to match
              annotations assembly using liftOver tool

  -x N        apply filtering of increased stringency when mapping SNP
              genomic positions to genes, N=0..3:
                0 - include all knownGene transcripts, no filtering
                1 - include only knownCanonical transcripts (default)
                2 - include only knownCanonical transcripts which
                    are also annotated as part of CCDS
                3 - include only knownCanonical transcripts with
                    exact match to CCDS

  -u          attempt to map translated CDS to a sequence in the local
              UniProtKB database, default is to utilize UCSC mappings

  -U          as above but utilize (a much slower) BLAT search instead
              of the exact matches; successful match requires 97%
              sequence identity and no more than 5 alignement gaps

  -b N        only include single transcript per each unique UniProtKB
              entry with the best matching CDS/protein sequences
              (requires -u|-U option), N=0..3:
                0 - accept all redundant UniProtKB matches (default)
                1 - all unique Swiss-Prot+TrEMBL+VarSplic entries
                2 - Swiss-Prot+TrEMBL entries (no VarSplic isoforms)
                3 - Swiss-Prot "canonical" entries only

  -s file     save to a FASTA-formatted file protein sequences for all
              CDS's to which at least a single missense SNP was mapped
              but no UniProtKB entries with matching sequences found

  -r N/M      partition input into M equally sized stripes and process
              stripe N

  -v N        verbosity level; default N=1, set N=0 for quiet operation

=head1 DESCRIPTION

Input file should be in tab-delimited text format with one genomic SNP
specification per each line, in the following format:

chrom:pos a1/a2[/a3[/a4]] [a1/a2 ...]

   * chrom  - chromosome identifier, e.g. chr10
   * pos    - SNP chromosome position, 1-based
   * a1/a2  - reference/observed allelic variants (nucleotides),
              specified on the plus strand of the assembly

Example:

   chr1:1267483   G/A
   chr1:1158631   A/C/G/T
   chr2:167262274 C/T
   chr4:264904    G/A
   chr7:122261636 G/T
   chr16:53698869 T/C

=head1 AUTHORS

Ivan Adzhubey

=head1 SUBVERSION

 $LastChangedDate: 2012-03-04 19:09:14 -0500 (Sun, 04 Mar 2012) $
 $LastChangedRevision: 395 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Getopt::Std;
use Pod::Usage;
use DBI;

my $PROGRESSBAR;

BEGIN {
  $Getopt::Std::STANDARD_HELP_VERSION = 1;
  use FindBin qw($RealBin);
  BEGIN {
    unless ($ENV{'PPH'}) {
      my $binpath = $RealBin;
      $binpath =~ s|/[^/]+$||;
      $ENV{'PPH'} = $binpath;
    }
  }
  use lib "$ENV{PPH}/perl";
  eval { require Term::ProgressBar; };
  $PROGRESSBAR = $@ ? 0 : 1;
}

use PPH::Config;
use PPH::Seq;

my $LD_PCRE='/usr/lib64/sqlite3/pcre.so';

my @TRANS = (
  # Standard genetic code
  {
    TTT => 'F',   TCT => 'S',   TAT => 'Y',   TGT => 'C',
    TTC => 'F',   TCC => 'S',   TAC => 'Y',   TGC => 'C',
    TTA => 'L',   TCA => 'S',   TAA => '*',   TGA => '*',
    TTG => 'L',   TCG => 'S',   TAG => '*',   TGG => 'W',

    CTT => 'L',   CCT => 'P',   CAT => 'H',   CGT => 'R',
    CTC => 'L',   CCC => 'P',   CAC => 'H',   CGC => 'R',
    CTA => 'L',   CCA => 'P',   CAA => 'Q',   CGA => 'R',
    CTG => 'L',   CCG => 'P',   CAG => 'Q',   CGG => 'R',

    ATT => 'I',   ACT => 'T',   AAT => 'N',   AGT => 'S',
    ATC => 'I',   ACC => 'T',   AAC => 'N',   AGC => 'S',
    ATA => 'I',   ACA => 'T',   AAA => 'K',   AGA => 'R',
    ATG => 'M',   ACG => 'T',   AAG => 'K',   AGG => 'R',

    GTT => 'V',   GCT => 'A',   GAT => 'D',   GGT => 'G',
    GTC => 'V',   GCC => 'A',   GAC => 'D',   GGC => 'G',
    GTA => 'V',   GCA => 'A',   GAA => 'E',   GGA => 'G',
    GTG => 'V',   GCG => 'A',   GAG => 'E',   GGG => 'G',
  },
  # Vertebrate mitochondrial code
  #
  # Differences from the Standard code:
  #         Code 2        Standard
  #  AGA    Ter  *          Arg  R
  #  AGG    Ter  *          Arg  R
  #  ATA    Met  M          Ile  I
  #  TGA    Trp  W          Ter  *
  {
    TTT => 'F',   TCT => 'S',   TAT => 'Y',   TGT => 'C',
    TTC => 'F',   TCC => 'S',   TAC => 'Y',   TGC => 'C',
    TTA => 'L',   TCA => 'S',   TAA => '*',   TGA => 'W',
    TTG => 'L',   TCG => 'S',   TAG => '*',   TGG => 'W',

    CTT => 'L',   CCT => 'P',   CAT => 'H',   CGT => 'R',
    CTC => 'L',   CCC => 'P',   CAC => 'H',   CGC => 'R',
    CTA => 'L',   CCA => 'P',   CAA => 'Q',   CGA => 'R',
    CTG => 'L',   CCG => 'P',   CAG => 'Q',   CGG => 'R',

    ATT => 'I',   ACT => 'T',   AAT => 'N',   AGT => 'S',
    ATC => 'I',   ACC => 'T',   AAC => 'N',   AGC => 'S',
    ATA => 'M',   ACA => 'T',   AAA => 'K',   AGA => '*',
    ATG => 'M',   ACG => 'T',   AAG => 'K',   AGG => '*',

    GTT => 'V',   GCT => 'A',   GAT => 'D',   GGT => 'G',
    GTC => 'V',   GCC => 'A',   GAC => 'D',   GGC => 'G',
    GTA => 'V',   GCA => 'A',   GAA => 'E',   GGA => 'G',
    GTG => 'V',   GCG => 'A',   GAG => 'E',   GGG => 'G',
  },
);

# Degenerate positions in codons:
#
# Nei M & Kumar S (2000) Molecular Evolution and Phylogenetics. Oxford University Press, New York.
# http://www.megasoftware.net/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_kumar_comeron_method.htm
# Note: in the following tables, '3' denotes a complex twofold site.
#
my @DGEN = (
  # Standard genetic code
  {
    TTT => [0,0,2],   TCT => [0,0,4],   TAT => [0,0,2],   TGT => [0,0,2],
    TTC => [0,0,2],   TCC => [0,0,4],   TAC => [0,0,2],   TGC => [0,0,2],
    TTA => [2,0,2],   TCA => [0,0,4],   TAA => [0,2,2],   TGA => [0,2,0],
    TTG => [2,0,2],   TCG => [0,0,4],   TAG => [0,0,2],   TGG => [0,0,0],

    CTT => [0,0,4],   CCT => [0,0,4],   CAT => [0,0,2],   CGT => [0,0,4],
    CTC => [0,0,4],   CCC => [0,0,4],   CAC => [0,0,2],   CGC => [0,0,4],
    CTA => [2,0,4],   CCA => [0,0,4],   CAA => [0,0,2],   CGA => [3,0,4],
    CTG => [2,0,4],   CCG => [0,0,4],   CAG => [0,0,2],   CGG => [3,0,4],

    ATT => [0,0,3],   ACT => [0,0,4],   AAT => [0,0,2],   AGT => [0,0,2],
    ATC => [0,0,3],   ACC => [0,0,4],   AAC => [0,0,2],   AGC => [0,0,2],
    ATA => [0,0,3],   ACA => [0,0,4],   AAA => [0,0,2],   AGA => [3,0,2],
    ATG => [0,0,0],   ACG => [0,0,4],   AAG => [0,0,2],   AGG => [3,0,2],

    GTT => [0,0,4],   GCT => [0,0,4],   GAT => [0,0,2],   GGT => [0,0,4],
    GTC => [0,0,4],   GCC => [0,0,4],   GAC => [0,0,2],   GGC => [0,0,4],
    GTA => [0,0,4],   GCA => [0,0,4],   GAA => [0,0,2],   GGA => [0,0,4],
    GTG => [0,0,4],   GCG => [0,0,4],   GAG => [0,0,2],   GGG => [0,0,4],
  },
  # Vertebrate mitochondrial code
  {
    TTT => [0,0,2],   TCT => [0,0,4],   TAT => [0,0,2],   TGT => [0,0,2],
    TTC => [0,0,2],   TCC => [0,0,4],   TAC => [0,0,2],   TGC => [0,0,2],
    TTA => [2,0,2],   TCA => [0,0,4],   TAA => [0,0,2],   TGA => [0,0,2],
    TTG => [2,0,2],   TCG => [0,0,4],   TAG => [0,0,2],   TGG => [0,0,2],

    CTT => [0,0,4],   CCT => [0,0,4],   CAT => [0,0,2],   CGT => [0,0,4],
    CTC => [0,0,4],   CCC => [0,0,4],   CAC => [0,0,2],   CGC => [0,0,4],
    CTA => [2,0,4],   CCA => [0,0,4],   CAA => [0,0,2],   CGA => [0,0,4],
    CTG => [2,0,4],   CCG => [0,0,4],   CAG => [0,0,2],   CGG => [0,0,4],

    ATT => [0,0,2],   ACT => [0,0,4],   AAT => [0,0,2],   AGT => [0,0,2],
    ATC => [0,0,2],   ACC => [0,0,4],   AAC => [0,0,2],   AGC => [0,0,2],
    ATA => [0,0,2],   ACA => [0,0,4],   AAA => [0,0,2],   AGA => [0,0,2],
    ATG => [0,0,2],   ACG => [0,0,4],   AAG => [0,0,2],   AGG => [0,0,2],

    GTT => [0,0,4],   GCT => [0,0,4],   GAT => [0,0,2],   GGT => [0,0,4],
    GTC => [0,0,4],   GCC => [0,0,4],   GAC => [0,0,2],   GGC => [0,0,4],
    GTA => [0,0,4],   GCA => [0,0,4],   GAA => [0,0,2],   GGA => [0,0,4],
    GTG => [0,0,4],   GCG => [0,0,4],   GAG => [0,0,2],   GGG => [0,0,4],
  },
);

# 0 - transition, 1 - transversion
my %TRNV = (
  'A' => { 'A' => undef, 'C' => 1, 'G' => 0, 'T' => 1 },
  'C' => { 'A' => 1, 'C' => undef,'G' => 1, 'T' => 0  },
  'G' => { 'A' => 0, 'C' => 1, 'G' => undef, 'T' => 1 },
  'T' => { 'A' => 1, 'C' => 0, 'G' => 1, 'T' => undef }
);

# Total 37 columns
my @header = qw(
  snp_pos
  str
  gene
  transcript
  ccid
  ccds
  cciden
  refa
  type
  ntpos
  nt1
  nt2
  flanks
  trv
  cpg
  jxdon
  jxacc
  exon
  cexon
  jxc
  dgn
  cdnpos
  frame
  cdn1
  cdn2
  aa1
  aa2
  aapos
  spmap
  spacc
  spname
  refs_acc
  dbrsid
  dbobsrvd
  dbavHet
  dbavHetSE
  dbRmPaPt
);

my @format = qw(
  %-16s
  %3s
  %12s
  %12s
  %5s
  %12s
  %6s
  %4s
  %12s
  %8s
  %3s
  %3s
  %6s
  %3s
  %3s
  %7s
  %7s
  %7s
  %7s
  %3s
  %3s
  %8s
  %5s
  %4s
  %4s
  %3s
  %3s
  %8s
  %5s
  %10s
  %12s
  %12s
  %10s
  %8s
  %10s
  %10s
  %8s
);

my $formatstr = join("\t", @format);
my $headerstr = sprintf('#'.$formatstr, @header) . "\tcomments\n";

$|++;

my %opts;
pod2usage('ERROR: No arguments specified') unless @ARGV;
getopts('aAchmnuUb:g:i:r:s:v:x:y:', \%opts) or pod2usage('ERROR: Wrong argument(s)');
pod2usage(-verbose=>2) if $opts{'h'};
pod2usage('ERROR: No input file specified') unless @ARGV;
my $infile = shift;

# Process options
my $VERBOSE = (exists $opts{'v'} && length $opts{'v'}) ? $opts{'v'} : 1;
my $target_assembly = $opts{'g'} || $CONFIG{'GENESET'} || 'hg19';
die "ERROR: Unsupported assembly version specified: $target_assembly\n"
  unless $target_assembly eq 'hg18' || $target_assembly eq 'hg19';
my $query_assembly  = $opts{'i'} || $target_assembly;
die "ERROR: Unsupported assembly version specified: $query_assembly\n"
  unless $query_assembly eq 'hg18' || $query_assembly eq 'hg19';
my $DATAPATH = "$CONFIG{GOLDENPATH}/$target_assembly/genes";
my $NUCFILE  = "$DATAPATH/knownGeneNuc.2bit";
my $KGDB     = "$DATAPATH/knownGene.sqlite";
my $MAPDB    = "$DATAPATH/kgToUp.sqlite";
my $CDSDB    = "$DATAPATH/knownGeneAA";
my $SNPDB    = "$DATAPATH/snp.sqlite";
my $coding_only   = $opts{'c'} || 0;
my $missense_only = $opts{'m'} || 0;
my $skip_noref    = $opts{'n'} || 0;
my $map2uniprot   = $opts{'u'} || 0;
$map2uniprot = 2 if $opts{'U'}; # super-meticious UniProt mapping using BLAT
my $stringency    = exists $opts{'x'} ? $opts{'x'} : 1;
die "ERROR: Illegal option argument: -x $stringency\n" unless $stringency =~ /^[0-3]$/;
# Only include single best unique UniProtKB entry hit per each transcript:
#  -b 1 - Swiss-Prot + TrEMBL + VarSplic unique entries
#  -b 2 - Swiss-Prot + TrEMBL entries (no VarSplic isoforms)
#  -b 3 - Swiss-Prot "canonical" entries only
my $besthit = (exists $opts{'b'} && length $opts{'b'}) ? $opts{'b'} : 0;
die "ERROR: Illegal option argument: -b $besthit\n" unless $besthit =~ /^[0-3]$/;
die "ERROR: -b option requires -u or -U also specified\n" if $besthit && ! $map2uniprot;
my $pph2output    = $opts{'y'} || 0;
my $seqfile = $opts{'s'} || '';
$CONFIG{'GRID_ARRAY'} = 1 if $opts{'a'};
$CONFIG{'GRID_ARRAY'} = 0 if $opts{'A'};
my ($pph2formatstr, $pph2headerstr);
if ($pph2output) {
  open(PPH2OUT, ">$pph2output") or die "ERROR: Can't create file: $pph2output\n";
  $pph2formatstr = join("\t", qw{ %-12s %6s %3s %3s %12s %6s %5s %3s %3s %4s %3s %3s %6s %6s });
  $pph2headerstr = sprintf('#'.$pph2formatstr, qw{ acc pos aa1 aa2 txname cdnpos frame nt1 nt2 5'3' trv cpg jxdon jxacc}) . "\tcomments\n";
}
$coding_only = $missense_only if $missense_only;

my $selectKg;
if ($besthit) {
  $selectKg = q{SELECT knownGene.*, macc, mcan FROM knownGene JOIN map ON name = qacc};
} else {
  $selectKg = q{SELECT * FROM knownGene};
}
$selectKg .= q{ WHERE exonCountCds > 0 AND chrom = ? AND txStart <= ? AND txEnd >= ?};
$selectKg .= q{ AND clusterId IS NOT NULL}  if $stringency >= 1;
$selectKg .= q{ AND ccdsId IS NOT NULL}     if $stringency >= 2;
$selectKg .= q{ AND cdsSimilarity = 1}      if $stringency >= 3;
$selectKg .= q{ AND}                        if $besthit > 1;
if      ($besthit == 0) {
  $selectKg .= q{ ORDER BY name};
} elsif ($besthit == 1) {
  $selectKg .= q{ ORDER BY mcan DESC, qfrac DESC, hits ASC, ident DESC, name};
} elsif ($besthit == 2) {
  $selectKg .= q{ NOT macc REGEXP '-\d+$' ORDER BY mcan DESC, qfrac DESC, hits ASC, ident DESC, name};
} elsif ($besthit == 3) {
  $selectKg .= q{ mcan = 1 ORDER BY qfrac DESC, hits ASC, ident DESC, name};
}

my $selectExChr = q{SELECT exStart, exEnd FROM kgExChr WHERE id = ? ORDER BY exStart ASC};
my $selectExCds = q{SELECT exFrom, exTo FROM kgExCds WHERE id = ? ORDER BY exFrom ASC};

my $selectMap = q{SELECT * FROM map JOIN hits USING(id) WHERE qacc = ?};
# Map only to the already preselected best-hit protein
$selectMap   .= q{ AND macc = ?} if $besthit;
$selectMap   .= q{ AND ? BETWEEN qfrom AND qto ORDER BY};
$selectMap   .= q{ mcan DESC,} unless $besthit;
$selectMap   .= q{ qfrac DESC, hits ASC, ident DESC LIMIT 1};

# Order in such a way that puts SNPs without valid orthologous nucleotides at the bottom
my $selectSnp = q{SELECT s.*, chimpAllele, orangAllele, macaqueAllele
FROM snp AS s LEFT JOIN snpOrtho AS o
USING(chrom,chromStart,strand,name)
WHERE chrom = ? AND chromStart = ?
ORDER BY chimpAllele IS NULL, chimpAllele = '?'};

my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $dbhKg  = DBI->connect("dbi:SQLite:dbname=$KGDB", '', '', $dbargs);
my $dbhMap = DBI->connect("dbi:SQLite:dbname=$MAPDB", '', '', $dbargs) if $map2uniprot;
my $dbhSnp = DBI->connect("dbi:SQLite:dbname=$SNPDB", '', '', $dbargs);

# attach kgToUp database to the same connection with knownGene
$dbhKg->do("ATTACH '$MAPDB' AS kgup") if $besthit;
# enable and load pcre extension module used by REGEXP operator to filter out isoforms
if ($besthit == 2 && $LD_PCRE) {
  my $dbh = DBI->connect("dbi:SQLite:dbname=", '', '', $dbargs);
  $dbh->sqlite_enable_load_extension(1);
  $dbh->do(qq{SELECT load_extension('$LD_PCRE')});
  $dbh->disconnect;
}

my $sthKg    = $dbhKg->prepare($selectKg);
my $sthExChr = $dbhKg->prepare($selectExChr);
my $sthExCds = $dbhKg->prepare($selectExCds);
my $sthMap   = $dbhMap->prepare($selectMap) if $map2uniprot;
my $sthSnp   = $dbhSnp->prepare($selectSnp);

# Check if we are running under SGE or LSF and if this is an array job
my ($grid_mode, $jobid, $taskid, $maxtask, $grid_msg);
if (exists $ENV{'SGE_ROOT'} && exists $ENV{'JOB_ID'}) {       # SGE
  $grid_mode = 'SGE';
  $jobid  = $ENV{'JOB_ID'} ne 'undefined' ? $ENV{'JOB_ID'} : 0;
  $taskid = $CONFIG{'GRID_ARRAY'} && exists $ENV{'SGE_TASK_ID'} && $ENV{'SGE_TASK_ID'} ne 'undefined' ? $ENV{'SGE_TASK_ID'} : 0;
} elsif ($ENV{'LSB_QUEUE'} && exists $ENV{'LSB_JOBID'}) {     # LSF
  $grid_mode = 'LSF';
  $jobid  = $ENV{'LSB_JOBID'};
  $taskid = $CONFIG{'GRID_ARRAY'} && exists $ENV{'LSB_JOBINDEX'} ? $ENV{'LSB_JOBINDEX'} : 0;
}
my $array_mode = $jobid && $taskid ? 1 : 0;

if ($grid_mode) {
  # Adjust PATH for a grid job
  if ($CONFIG{'GRID_PATH'}) {
    my $ARCH = `uname -m`; chomp $ARCH;
    $ENV{'PATH'} = eval "qq{$CONFIG{GRID_PATH}}";
  }
  $grid_msg = "I am $grid_mode job $jobid \@ $ENV{HOSTNAME}\n";
}

if ($array_mode) {
  # GGI_STRIPES overwrites any array size specified initially
  $maxtask = $ENV{'GGI_STRIPES'} || $ENV{'SGE_TASK_LAST'} || $ENV{'LSB_JOBINDEX_END'};
  $grid_msg = "I am $grid_mode job $jobid task $taskid of $maxtask \@ $ENV{HOSTNAME}\n";
}

my ($stripe_no, $stripes);
if ($opts{'r'}) {
  ($stripe_no, $stripes) = $opts{'r'} =~ m|^(\d+)/(\d+)$|;
  die "ERROR: Illegal option's argument: -r $opts{'r'}\n" unless $stripe_no && $stripes;
}

if ($seqfile) {
  warn "WARNING: File exists, overwriting: $seqfile\n" if -e $seqfile && -s $seqfile;
  open(SEQFILE, ">$seqfile") or die "Can't create file: $seqfile\n";
}

# Debugging
warn $grid_msg if defined $grid_msg;

if ($array_mode) {
  warn "WARNING: Explicit input stripe specification via '-r' option ignored in $grid_mode array mode\n" if $opts{'r'};
  # Additional sanity checks specific to PPHWeb2 service:
  # In the web service pipeline mode (under apache2 user), a special environment
  # variable has to be set since SGE_TASK_LAST cannot be altered for a pending job
  if ($grid_mode eq 'SGE' && $ENV{'USER'} eq 'apache2') {
    die "ERROR: Missing GGI_STRIPES environment variable while in SGE array mode\n"
      unless exists $ENV{'GGI_STRIPES'} && length $ENV{'GGI_STRIPES'};
    die "ERROR: Illegal GGI_STRIPES value while in SGE array mode: $ENV{GGI_STRIPES}\n"
      unless $ENV{'GGI_STRIPES'} > 0;
  }
  ($stripe_no, $stripes) = ($taskid, $maxtask);
}

warn "Initializing ...\n" if $VERBOSE;

# Count how many input lines we have and how many of these are valid query lines,
# e.g., sans comments and empty lines. Used by progress bar and with the stripes option.
my $input_line_count  = 0;
my $input_query_count = 0;
if ($VERBOSE || $stripes) {
  open(FIN, $infile) or die "ERROR: Can't open input file: $infile\n";
  while (<FIN>) {
    $input_line_count++;
    # Skip comments and empty lines
    next if /^#/ || /^\s*$/;
    $input_query_count++;
  }
  close(FIN);
}

# Calculate stripe range
my ($range_start, $range_end) = stripe_range($stripe_no, $stripes, $input_query_count) if $stripes;

my %totals;
my $limit = 0;
my $input_line_no  = 0;
my $input_query_no = 0;
my $range_line_no  = 0;
my %genes;
my $input_range_count;
if      ($stripes) {
  if ($range_start == 0 && $range_end == 0) {
    $input_range_count = 0;
  } else {
    $input_range_count = $range_end - $range_start + 1
  }
  die "ERROR: Input out of range (stripe $stripe_no / $stripes with $input_query_count query lines in input)\n"
    unless $input_range_count;
} elsif ($VERBOSE) {    # $input_query_count is only set in either -r/array or verbose modes
  $input_range_count = $input_query_count;
  die "ERROR: Nothing to do ($input_query_count query lines in input)\n"
    unless $input_range_count;
}
my $progress = Term::ProgressBar->new({count=>$input_range_count}) if $PROGRESSBAR && $VERBOSE;  # set up a progress bar
my $next_update        = 0;
my $header_printed     = 0;
my $pph2header_printed = 0;
open(FIN, $infile) or die "ERROR: Can't open input file: $infile\n";
LINE: while (<FIN>) {
  $input_line_no++;   # unfiltered input line number, used in error/warning messages
  # Skip comments and empty lines
  next if /^#/ || /^\s*$/;
  $input_query_no++;  # input query line number
  if ($stripes) {
    next LINE if $input_query_no < $range_start;
    last LINE if $input_query_no > $range_end;
  }
  $totals{'input'}++; # total query input lines to process (i.e., within range but without comments)
  $range_line_no++;   # count, same as above
  $next_update = $progress->update($range_line_no) if $PROGRESSBAR && $VERBOSE && $input_line_no >= $next_update;
  chomp;
  # Process optional user comments as separated by a '#' character
  my $comments = $1 if s/\s *#\s*(.*)$//;
  my ($chrpos, @a) = split;
  my ($chrom, $snppos);
  if ($chrpos =~ /^chr(\d{1,2}|[MXY]):(\d+)/i) {
    $chrom  = 'chr' . $1;
    $snppos = $2;
  } else {
    warn "WARNING: Illegal genomic coordinate specification skipped at input line ($input_line_no): $chrpos\n";
    next;
  }
  my @o_alleles;   # these are original alleles specified in input,
  my %o_alleles;   # they are always on plus strand
  for (my $i=0; $i<=$#a; $i++) {
    $a[$i] =~ s/\s+//g;
    warn("WARNING: Illegal allele specification ($a[$i]) skipped at input line ($input_line_no)\n"), next
      unless $a[$i] =~ m|^[ACGTacgt]/[ACGTacgt](?:[/,][ACGTacgt])*$|;
    foreach my $nuc (split m|[/,]|, $a[$i]) {
      $nuc = uc $nuc;
      push @o_alleles, $nuc unless $o_alleles{$nuc};
      $o_alleles{$nuc} = 1;
    }
  }
  warn("WARNING: No valid alleles found at input line ($input_line_no)\n"), next LINE unless @o_alleles;
  last LINE if $limit and $input_line_no > $limit;

  my $snpstart = $snppos   - 1; # convert to zero base coordinate
  my $snpend   = $snpstart + 1; # dummy open end coordinate

  # liftOver input coordinates to target assembly if query uses different assembly version
  if ($query_assembly ne $target_assembly) {
    my $fhq = File::Temp->new();
    my $fnq = $fhq->filename;     # query in BED format
    my $fht = File::Temp->new();
    my $fnt = $fht->filename;     # target in BED format
    my $fhu = File::Temp->new();
    my $fnu = $fhu->filename;     # unMapped
    print($fhq "$chrom\t$snpstart\t$snpend\n");
    my $lofile =
      "$CONFIG{GOLDENPATH}/$query_assembly/liftOver/" .
      $query_assembly . 'To' . ucfirst($target_assembly) .'.over.chain';
    die "ERROR: liftOver file not found: $lofile\n" unless -e $lofile;
    my $rc = system('liftOver', $fnq, $lofile, $fnt, $fnu) == 0 or
      warn("WARNING: liftOver failed at input line ($input_line_no): $!\n"), next LINE;
    my ($mapped, $unmapped);
    { local $/; $mapped   = <$fht>; }
    { local $/; $unmapped = <$fhu>; }
    chomp($mapped, $unmapped);
    if (length $unmapped) {
      $unmapped =~ s/^#//s;
      $unmapped =~ s/\n/: /mg;
      $unmapped =~ s/\t+/ /mg;
      warn "WARNING: liftOver error at input line ($input_line_no): $unmapped\n";
      next LINE;
    } else {
      warn("WARNING: liftOver output empty at input line ($input_line_no)\n"), next LINE
        unless length $mapped;
      ($chrom, $snpstart, $snpend) = split ' ', $mapped;
    }
  }

  $totals{'processed'}++; # total query input lines processed (after validating & without comments)

  # Execute query vs knownGene table of gene transcripts
  $sthKg->execute($chrom, $snpstart, $snpend);

  # Standard genetics code index is 0, vertebrate mitochondrial code is 1
  my $gencode = $chrom eq 'chrM' ? 1 : 0;

  my (%qacclist, %macclist);
  GENE: while (my $knownGene = $sthKg->fetchrow_hashref) {

    if ($besthit) {
      # skip all but the very first/best unique transcript to unique protein hit (-b1)
      next GENE if $qacclist{$knownGene->{name}} || $macclist{$knownGene->{macc}};
      # register this transcript/protein hit
      $qacclist{$knownGene->{name}} = $macclist{$knownGene->{macc}} = 1;
    }

    my $nucseq;
    if (exists $genes{$knownGene->{name}}{txseq}) {
      $nucseq = $genes{$knownGene->{name}}{txseq};
    } else {
      $nucseq = PPH::Seq::fetch_txNuc($knownGene->{name}, $NUCFILE);
      $genes{$knownGene->{name}}{txseq} = $nucseq;
    }
    next GENE unless defined $nucseq;

    my $minusflag = $knownGene->{strand} eq '-' ? 1 : 0;

    # Extract chromosome nt at the variation position and at its 5' and 3' flanks (for CpG annotation)
    # and complement working copies of input alleles (array and hash) for minus strand transcripts.
    # Note1: all three nt's are on the coding strand, with flanks enumerated in the direction of transcription (5'3')
    # Note2: corresponding flanking nt will be undef for a SNP at the first / last positions of the transcript
    my @alleles = @o_alleles;   # these are working copies of the original input alleles, they will be
    my %alleles = %o_alleles;   # complemented for positions mapped to transcripts on minus strand
    # SNP position in the transcript nucleotide sequence, zero base
    my $nucpos  = $minusflag ? $knownGene->{txEnd} - $snpend : $snpstart - $knownGene->{txStart};
    # Transcript nt at the SNP position
    my $tx_nt   = substr($nucseq, $nucpos, 1);
    # Reference assembly nt at the SNP position, will be complemented later for transcripts on minus strand
    my $rf_nt   = $tx_nt;
    unless (defined $tx_nt) {
      warn "ERROR: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Failed to extract reference nucleotide\n";
      next GENE;
    }
    unless ($tx_nt =~ /^[ACTG]$/) {
      warn "ERROR: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Illegal reference nucleotide ($tx_nt) in transcript sequence\n";
      next GENE;
    }
    my $nuclen = length $nucseq;
    my $tx_nt5 = $nucpos > 0           ? substr($nucseq, $nucpos - 1, 1) : '?';
    my $tx_nt3 = $nucpos < $nuclen - 1 ? substr($nucseq, $nucpos + 1, 1) : '?';
    if ($minusflag) {
      # Complement back to plus (reference) strand of assembly
      $rf_nt   =~ tr/ACGT/TGCA/;
      # Complement input alleles to minus (transcript coding) strand
      %alleles = ();
      for (my $i=0; $i<@alleles; $i++) {
        $alleles[$i] =~ tr/ACGT/TGCA/;
        $alleles{$alleles[$i]} = 1;
      }
    }
    # Check for the situation when reference nt is different from any of the variant alleles submitted by the user
    unless ($o_alleles{$rf_nt}) {
      warn("WARNING: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) None of the input alleles (".join('/', @o_alleles).") matches reference allele ($rf_nt)\n");
      $totals{'processed'}--, next LINE if $skip_noref;
    }

    # Scan exons to obtain SNP exon-based annotations (SNP functional class,
    # SNP position in CDS, exon/intron junction distances, etc)
    my $snptype;
    my ($txexno, $cdsexno);
    my ($excount, $cdscount) = (0, 0);
    my ($jxdon, $jxacc);  # 5' donor (>GT), 3' acceptor (AG<)
    my ($dstart, $dend);
    my $cdslen = 0;       # CDS length (number of nucleotides)
    my $cdspos;           # SNP position in the CDS nucleotide sequence, zero base
    $sthExChr->execute($knownGene->{id});
    while (my $exon = $sthExChr->fetchrow_hashref) {

      my $exstart = $exon->{exStart};
      my $exend   = $exon->{exEnd};

      $excount++;

      # Annotate SNP position
      unless ($snptype) {
        # SNP position is inside this exon
        if ($snpstart >= $exstart && $snpend <= $exend) {
          $txexno = $excount;
          # SNP position is in 5'-UTR
          if      ($snpstart < $knownGene->{cdsStart}) {
            $snptype = $minusflag ? 'utr-3' : 'utr-5';
          # SNP position is in 3'-UTR
          } elsif ($snpend   > $knownGene->{cdsEnd})    {
            $snptype = $minusflag ?  'utr-5' : 'utr-3';
          # SNP position is in the coding exon sequence
          } else {
            $snptype = 'exon';
            $cdsexno = $cdscount + 1;
            $cdspos  = $cdslen + $snpstart -
                       ($exstart < $knownGene->{cdsStart} ? $knownGene->{cdsStart} : $exstart);
          }
        }
      }

      # Locate nearest exon/intron boundaries
      if ($knownGene->{exonCount} > 1) {
        # Skip exon start if it is at the transcript start -- no junction
        unless ($exstart == $knownGene->{txStart}) {
          my $ds  = $exstart - $snpstart - 1;
          $dstart = $ds if !defined $dstart || abs($ds) < abs($dstart);
        }
        # Skip exon end if it is at the transcript end -- no junction
        unless ($exend == $knownGene->{txEnd}) {
          my $de  = $exend - $snpend + 1;
          $dend   = $de if !defined $dend || abs($de) < abs($dend);
        }
      }

      # Skip completely non-coding exons
      next if $exstart >= $knownGene->{cdsEnd} || $exend <= $knownGene->{cdsStart};

      $cdscount++;

      $cdslen += ($exend   > $knownGene->{cdsEnd}   ? $knownGene->{cdsEnd}   : $exend)  -
                 ($exstart < $knownGene->{cdsStart} ? $knownGene->{cdsStart} : $exstart);

    }

    if ($cdslen % 3) {
      warn "ERROR: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Incomplete CDS (len=$cdslen)\n";
      next GENE;
    }

    # Revert SNP-containing exon numbers to match the transcript direction
    if ($minusflag) {
      $txexno  = $knownGene->{exonCount}    - $txexno  + 1 if defined $txexno;
      $cdsexno = $knownGene->{exonCountCds} - $cdsexno + 1 if defined $cdsexno;
    }

    $snptype = 'intron' unless defined $snptype;
    $txexno  = '?'      unless defined $txexno;
    $cdsexno = '?'      unless defined $cdsexno;

    if (defined $dstart && defined $dend) {
      # dstart and dend are swapped and reversed for minus strand transcripts
      ($dstart, $dend) = (-$dend, -$dstart) if $minusflag;
      $jxdon = sprintf('%+d', $dend);
      $jxacc = sprintf('%+d', $dstart);
    }

    # Annotate all coding SNP alleles
    if ($snptype eq 'exon') {

      $cdspos = $cdslen - $cdspos - 1 if $minusflag;
      my $cdnpos = int($cdspos/3);  # CDS codon number where the SNP is located, zero base
      my $frame  = $cdspos % 3;     # SNP position within the codon, [0..2]
      my $cdnnuc = $cdnpos * 3;     # Position in CDS of the first nucleotide of codon where SNP is located, zero base

      # Maps nucleotide positions in the CDS sequence to full transcript sequence positions
      my @cds2gene;
      $sthExCds->execute($knownGene->{id});
      while (my $cdsExon = $sthExCds->fetchrow_hashref) {
        push @cds2gene, $cdsExon->{exFrom} .. $cdsExon->{exTo};
      }

      unless (defined $cds2gene[$cdnnuc] && defined $cds2gene[$cdnnuc+2]) {
        warn "ERROR: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Illegal transcript codon position ($cdnnuc) for exonic SNP\n";
        next GENE;
      }

      my $cdn1 =
        substr($nucseq, $cds2gene[$cdnnuc],   1) .
        substr($nucseq, $cds2gene[$cdnnuc+1], 1) .
        substr($nucseq, $cds2gene[$cdnnuc+2], 1);

      my $aa1 = $TRANS[$gencode]{$cdn1};

      # Skip ambiguous / unknown codons
      unless (defined $aa1) {
        warn "ERROR: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Ambiguous or illegal reference codon ($cdn1) in CDS at position ($cdspos)\n";
        next GENE;
      }

      # Check if the mutated codon is split by an exon/intron junction boundary
      #        d   a        (d - donor site; a - acceptor site)
      #   =====|---|ACG=     d << 0 / f(0),a(-1) f(1),a(-2) f(2),a(-3)
      #   ====A|---|CG==     f(0),d(+1) / f(1),a(-1) f(2),a(-2)
      #   ===AC|---|G===     f(0),d(+2) f(1),d(+1) / f(2),a(-1)
      #   ==ACG|---|====     f(0),d(+3) f(1)d(+2) f(2),d(+3) / a >> 0
      my $codon_split;
      if (defined $jxdon && defined $jxacc) {
        if      ($frame == 0) {
          $codon_split = 1 if $jxdon ==  1 || $jxdon ==  2;
        } elsif ($frame == 1) {
          $codon_split = 1 if $jxdon ==  1 || $jxacc == -1;
        } elsif ($frame == 2) {
          $codon_split = 1 if $jxacc == -1 || $jxacc == -2;
        }
      }

      # Obtain transcript->protein mapping if one of the -u or -U options was specified.
      # Note, this mapping cannot be cached because different CDS positions can map to
      # different UniProtKB sequences (such mappings are usually dubious but retained
      # in the map database in order to maximize coverage).
      my ($mapped, $aaoffset, $aapos, $mident);
      if ($map2uniprot) {
        if ($besthit) {
          $sthMap->execute($knownGene->{name}, $knownGene->{macc}, $cdnpos);
        } else {
          $sthMap->execute($knownGene->{name}, $cdnpos);
        }
        $mapped = $sthMap->fetchrow_hashref;
        $sthMap->finish;
        if (defined $mapped) {
          # Position offset between translated transcript CDS and matched UniProtKB protein sequences
          $aaoffset = $mapped->{mfrom} - $mapped->{qfrom};
          # Mapped codon / AA residue position in the UniProtKB protein sequence, zero base
          $aapos    = $cdnpos + $aaoffset;
          # Identity in overlap
          $mident   = $mapped->{ident};
        } else {
          $mident   = 0;
          warn "WARNING: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) No UniProtKB matches found for CDS codon position ($cdnpos)\n" if $VERBOSE > 1;
        }
      }

      ALLELE: for (my $i=0; $i<@alleles; $i++) {
        my $allele = $alleles[$i];
        next ALLELE if $allele eq $tx_nt;
        my $o_allele = $o_alleles[$i];
        my $cdn2 = $cdn1; substr($cdn2, $frame, 1, $allele);
        my $aa2  = $TRANS[$gencode]{$cdn2} if exists $TRANS[$gencode]{$cdn2};
        if      ($aa2 eq $aa1) {
          $snptype = 'coding-synon';
        } elsif ($aa2 eq '*' || $aa1 eq '*') {
          $snptype = 'nonsense';
        } else {
          $snptype = 'missense';
        }
        $totals{$snptype}++; $totals{'annotated'}++;
        next ALLELE if $missense_only && $snptype ne 'missense';
        my $flanks = $tx_nt5 . $tx_nt3;
        my $cpg    = CpG($tx_nt5, $tx_nt, $allele, $tx_nt3);
        my $dbsnp  = map2dbSNP($chrom, $snpstart, $tx_nt, $allele, $knownGene);
        # Save translated CDS protein sequence if -s option was specified
        if ($seqfile && $snptype eq 'missense' && !$genes{$knownGene->{name}}{saved}) {
          # Save the protein sequence if either no mapping was attempted and at the same time spID was
          # absent from kgXref table or mapping attempt has failed to find a matching UniProtKB sequence
          if ($map2uniprot && !defined $mapped || !$map2uniprot && !$knownGene->{spID}) {
            my ($head, $acc, $name, $desc, $seq) = PPH::Seq::fastacmd($knownGene->{name}, $CDSDB);
            print SEQFILE $head, "\n", $seq, "\n";
          }
          $genes{$knownGene->{name}}{saved} = 1;
        }
        # Print out complete SNP annotation
        my @values = map { (defined) ? $_ : '' } (
          "$chrom:$snpend",
          $knownGene->{strand},
          $knownGene->{geneSymbol},
          $knownGene->{name},
          $knownGene->{clusterId},
          $knownGene->{ccdsId},
          $knownGene->{cdsSimilarity},
          "$rf_nt/$o_allele",
          $snptype,
          defined $nucpos ? $nucpos + 1 : '',
          $tx_nt,
          $allele,
          $flanks,
          $TRNV{$tx_nt}{$allele},
          $cpg,
          $jxdon,
          $jxacc,
          "$txexno/$knownGene->{exonCount}",
          "$cdsexno/$knownGene->{exonCountCds}",
          $codon_split,
          $DGEN[$gencode]{$cdn1}[$frame],
          defined $cdnpos ? $cdnpos + 1 : '',
          $frame,
          $cdn1,
          $cdn2,
          $aa1,
          $aa2,
          $map2uniprot && defined $aapos ? $aapos + 1 : '',
          $map2uniprot ? $mident : '',
          $map2uniprot && defined $mapped ? $mapped->{macc} : $knownGene->{spID},
          $map2uniprot && defined $mapped ? $mapped->{mname} : $knownGene->{spDisplayID},
          $knownGene->{protAcc},
          $dbsnp->{name},
          $dbsnp->{observed},
          $dbsnp->{avHet},
          $dbsnp->{avHetSE},
          $dbsnp->{RmPaPt}
        );
        unless ($header_printed || $stripes && $stripe_no > 1) {
          print $headerstr; $header_printed = 1;
        }
        printf $formatstr, @values;
        print "\t# $comments" if defined $comments && length $comments;
        print "\n";
        # Print out query string for missense SNPs in PolyPhen-2 input format
        # to a filename specified via -y option.
        if ($pph2output && $snptype eq 'missense') {
          my $spID = $knownGene->{spID} if $knownGene->{spID};
          # Mapping to UniProtKB was attempted...
          if ($map2uniprot) {
            # ...and succeeded, replace kgXref spID with the UniProtKB one
            if (defined $mapped) {
              $spID = $mapped->{macc};
            # ...and failed, skip unless -s option was also specified,
            # in which case substitute UCSC transcript name for spID
            } else {
              if ($seqfile) { $spID = $knownGene->{name}; $aapos = $cdnpos; }
              else { $spID = ''; }
            }
          # No UniProtKB mapping, use original kgXref accession if available...
          } else {
            # ...otherwise skip unless -s option was also specified,
            # in which case substitute UCSC transcript name for spID
            if (!$spID && $seqfile) { $spID = $knownGene->{name}; }
            $aapos = $cdnpos if $spID;
          }
          if ($spID) {
            # Remove '/' character from original (input) allele spec since snp_id is
            # used in XML file names and slash is not allowed in filenames
            $o_allele =~ tr|/||d;
            my $snpid = "$chrom:$snpend|$rf_nt$o_allele|$knownGene->{name}$knownGene->{strand}|";
            $snpid   .= $knownGene->{geneSymbol} if defined $knownGene->{geneSymbol};
            $snpid   .= '|';
            $snpid   .= $knownGene->{protAcc} if defined $knownGene->{protAcc};
            my @values = map { (defined) ? $_ : '' } (
              $spID,
              $aapos + 1,
              $aa1,
              $aa2,
              $knownGene->{name},
              $cdnpos + 1,
              $frame,
              $tx_nt,
              $allele,
              $flanks,
              $TRNV{$tx_nt}{$allele},
              $cpg,
              $jxdon,
              $jxacc
            );
            unless ($pph2header_printed || $stripes && $stripe_no > 1) {
              print PPH2OUT $pph2headerstr; $pph2header_printed = 1;
            }
            printf PPH2OUT $pph2formatstr, @values;
            # PolyPhen-2 input format requires comments to be additionally delimited with the '#' character
            my $idcomm = defined $comments && length $comments ? "$snpid\t$comments" : $snpid;
            print PPH2OUT "\t# $idcomm\n";
          }
        }
      } # ALLELE

    # Short annotation for non-coding SNP alleles
    } else {
      ALLELE: for (my $i=0; $i<@alleles; $i++) {
        my $allele = $alleles[$i];
        next ALLELE if $allele eq $tx_nt;
        $totals{$snptype}++; $totals{'annotated'}++;
        # Skip reporting of all non-coding alleles if '-c' option was specified
        # but still count totals for them properly in the above line.
        next ALLELE if $coding_only;
        my $o_allele = $o_alleles[$i];
        my $flanks = $tx_nt5 . $tx_nt3;
        my $cpg   = CpG($tx_nt5, $tx_nt, $allele, $tx_nt3);
        my $dbsnp = map2dbSNP($chrom, $snpstart, $tx_nt, $allele, $knownGene);
        # Print out complete SNP annotation
        my @values = map { (defined) ? $_ : '' } (
          "$chrom:$snpend",
          $knownGene->{strand},
          $knownGene->{geneSymbol},
          $knownGene->{name},
          $knownGene->{clusterId},
          $knownGene->{ccdsId},
          $knownGene->{cdsSimilarity},
          "$rf_nt/$o_allele",
          $snptype,
          defined $nucpos ? $nucpos + 1 : '',
          $tx_nt,
          $allele,
          $flanks,
          $TRNV{$tx_nt}{$allele},
          $cpg,
          $jxdon,
          $jxacc,
          (undef) x 15,
          $dbsnp->{name},
          $dbsnp->{observed},
          $dbsnp->{avHet},
          $dbsnp->{avHetSE},
          $dbsnp->{RmPaPt}
        );
        unless ($header_printed || $stripes && $stripe_no > 1) {
          print $headerstr; $header_printed = 1;
        }
        printf $formatstr, @values;
        print "\t# $comments" if defined $comments && length $comments;
        print "\n";
      } # ALLELE
    }

  } # GENE

} # LINE
close(FIN);

$dbhSnp->disconnect;
$dbhMap->disconnect if $map2uniprot;
$dbhKg->disconnect;

close(SEQFILE) if $seqfile;
close(PPH2OUT) if $pph2output;

# Print out totals
my $skipped = ($totals{'input'} || 0) - ($totals{'processed'} || 0);
print "## Totals:\n";
printf "##   %-18s %8s\n", 'lines input', $totals{'input'} || 0;
printf "##   %-18s %8s\n", 'lines skipped', $skipped;
printf "##   %-18s %8s\n", 'alleles annotated', $totals{'annotated'} || 0;
delete $totals{'input'}; delete $totals{'processed'}; delete $totals{'annotated'};
for my $key (qw{ missense nonsense coding-synon intron utr-3 utr-5 }) {
  printf "##     %-16s %8s\n", $key, $totals{$key} || 0;
  delete $totals{$key};
}
# In case something non-standard was still left
for my $key (sort keys %totals) {
  printf "##     %-16s %8s\n", $key, $totals{$key} || 0;
}

warn "Finished.\n" if $VERBOSE;

#------------------------------------------------------------------------------

# Calculate [start,end] range of a stripe in dataset
# Note: all indices passed and returned are one-based
#
# Arguments:  stripe index, number of stripes, dataset size
# Returns:    data start index, data end index; both indices
#             returned are zero if no more data fall within
#             a stripe specified (e.g. number of stripes is
#             larger than the number of data items)
#
sub stripe_range {
  my ($stripe_no, $stripes, $size) = @_;

  die "ERROR: Stripe index too large: $stripe_no > $stripes\n" if $stripe_no > $stripes;

  my $p = int($size / $stripes);
  my $r = $size % $stripes;

  my $shift_start = ($stripe_no <= $r ? $stripe_no : $r + 1) - 1;
  my $shift_end   = ($stripe_no <= $r ? $stripe_no : $r);

  my $start = ($stripe_no - 1)*$p + 1 + $shift_start;
  my $end   = $stripe_no*$p + $shift_end;

  if ($start > $end) {
    return (0, 0);
  } else {
    return ($start, $end);
  }
}

sub CpG {
  my ($nt5, $allele1, $allele2, $nt3) = @_;
  my $t1  = $nt5 . $allele1 . $nt3;
  my $t2  = $nt5 . $allele2 . $nt3;
  warn("WARNING: Incorrect triplet(s): $t1>$t2\n"), return
    unless length($t1) == 3 && length($t2) == 3;
  my $cpg = 0;
  $cpg    = 1 if index($t1, 'CG') >= 0;
  $cpg   += 2 if index($t2, 'CG') >= 0;
  return $cpg;
}

# Note: snpXXXOrthoPt2Pa2Rm2.humanAllele contains a *reference* (plus strand) allele,
#       unlike the corresponding fields for all other species, for which alleles are
#       always relative to the human SNP strand. I assume this is a bug in UCSC data.
#       On the contrast, snpXXXOrthoPt2Pa2Rm2.humanObserved as well as snpXXX.observed
#       are listed relative to the current SNP strand, just like all other alleles.
sub map2dbSNP {
  my ($chrom, $snpstart, $tx_nt, $allele, $knownGene) = @_;
  my $snpend = $snpstart + 1;
  $sthSnp->execute($chrom, $snpstart);
  my $firstdbsnp;
  my @extrasnps;
  while (my $dbsnp = $sthSnp->fetchrow_hashref) {
    # Complement observed alleles if SNP orientation does not match gene orientation
    $dbsnp->{observed} =~ tr/ACGTacgt/TGCAtgca/ if $knownGene->{strand} ne $dbsnp->{strand};
    # Both reference and variant allele should match (in any order)
    unless ($dbsnp->{observed} =~ /$tx_nt/ && $dbsnp->{observed} =~ /$allele/) {
      warn "WARNING: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Allele ($tx_nt/$allele) does not match dbSNP ($dbsnp->{name}) observed ($dbsnp->{observed})\n" if $VERBOSE > 1;
      next;
    }
    # Skip all but first match
    if (defined $firstdbsnp) {
      push @extrasnps, $dbsnp->{name};
    } else {
      $firstdbsnp = $dbsnp;
    }
  }
  if (@extrasnps && $VERBOSE > 1) {
    warn "WARNING: ($chrom:$snpend $knownGene->{strand} $knownGene->{name}) Allele ($tx_nt/$allele) matches multiple dbSNPs: $firstdbsnp->{name}, " .
      join(', ', @extrasnps) . "\n";
  }
  if (defined $firstdbsnp) {
    $firstdbsnp->{avHet}   = (exists $firstdbsnp->{avHet} && defined $firstdbsnp->{avHet}) ?
                              sprintf('%g', $firstdbsnp->{avHet}) : '';
    $firstdbsnp->{avHetSE} = (exists $firstdbsnp->{avHetSE} && defined $firstdbsnp->{avHetSE}) ?
                              sprintf('%g', $firstdbsnp->{avHetSE}) : '';
    # If JOIN with the snpXXXOrthoPt2Pa2Rm2 table was successful then all alleles
    # are defined, even though any number of them could be ?, or unknown
    if (defined $firstdbsnp->{macaqueAllele}) {
      $firstdbsnp->{RmPaPt} =
        "$firstdbsnp->{macaqueAllele}>$firstdbsnp->{orangAllele}>$firstdbsnp->{chimpAllele}";
      # Complement all orthologs alleles if SNP orientation does not match gene orientation
      $firstdbsnp->{RmPaPt} =~ tr/ACGTacgt/TGCAtgca/ if $firstdbsnp->{strand} ne $knownGene->{strand};
    }
  }
  return $firstdbsnp;
}
