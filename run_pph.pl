#!/usr/bin/env perl

use warnings;
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=pod

=head1 NAME

run_pph.pl - run a set of PolyPhen-2 queries

=head1 SYNOPSIS

run_pph.pl [options] infile [>outfile]

where options are:

  -s seqfile  read query protein sequences from FASTA-formatted seqfile

  -b dbase    BLAST sequence database; default is: nrdb/uniref100

  -d dumpdir  directory for auxiliary output files; default is: scratch

  -r N/M      split input file into M parts with equal number of query
              lines each and process part N

  -v level    debugging output, verbosity level=[1..6]

=head1 DESCRIPTION

Input file should consist of tab-delimited lines with one amino acid
residue substitution per line, described by the following 4 mandatory
fields:

acc pos aa1 aa2

   * acc       - protein accession or entry name; sequences will be
                 retrieved from the local copy of UniProtKB database
                 or from an optional user-supplied FASTA file
   * pos       - substitution position in the protein sequence
   * aa1       - reference AA residue as found at the substitution
                 position in the query protein sequence
   * aa2       - substitution AA residue

Note: aa1 & aa2 are standard 1-letter residue codes

See sets/test.input file included for the example of input format.

Advanced options:

  -e eval   change E-value cutoff for BLAST searches; default is: 1e-3

  -p        use precomputed multiple alignments only

  -f        use full length protein sequences for multiple alignments;
            default is to use BLAST HSPs

  -t        create alignments with BLAST rather than MAFFT/LEON/Cluspack
            pipeline (use for debugging purposes only)

  -m        attempt to map proteins to genomic sequences to extract
            codons and flanking nucleotides

  -g hgXX   genome assembly version used when mapping proteins to genes;
            default is hg19 (GRCh37), also supported is hg18 (NCBI36)

  -x        save detailed SNP annotations to XML file in the current
            directory

  -c path   read configuration files from a directory path

=head1 AUTHORS

Ivan Adzhubey, Steffen Schmidt

=head1 SUBVERSION

 $LastChangedDate: 2021-05-23 20:32:04 -0400 (Sun, 23 May 2021) $
 $LastChangedRevision: 406 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

our $VERSION;

use Getopt::Std;
use Pod::Usage;
use File::Basename;
use File::Temp qw(tempdir);
use File::Path;
use Cwd qw(abs_path cwd);
use CGI;

BEGIN {
  $VERSION = '2.2.3r406';
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
}

# Need this before BEGIN block to import $DEBUG tested in it
use PPH::Config;

# Use BEGIN block here to ensure that command-line options are
# processed before Smart::Comments module is loaded since DEBUG
# flag can be altered from the command line.
my %opts;
BEGIN {
  # Process command-line options
  pod2usage('No arguments specified!') unless @ARGV;
  getopts('aAfhmptxzb:c:d:e:g:q:r:s:v:', \%opts) or pod2usage('Wrong arguments!');
  pod2usage(-verbose=>2) if $opts{'h'};
  pod2usage('No input file(s) specified!') unless $opts{'q'} || @ARGV;
  # Adjust default configuration if options were specified on command line
  PPH::Config->new(\%opts) if %opts;
  if ($DEBUG) {
    # Test if we have Smart::Comments installed
    eval { require Smart::Comments; Smart::Comments->import(-ENV) };
    die "Please install Smart::Comments module or switch off verbose/DEBUG mode\n" if $@;
  }
}

# All other PPH modules should follow PPH::Config for %CONFIG hash to be
# properly set before the modules that may attempt to utilize its contents.
use PPH;
use SNP;

$| = 1;

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
# Starting...
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

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
  $maxtask  = $ENV{'GGI_STRIPES'} || $ENV{'SGE_TASK_LAST'} || $ENV{'LSB_JOBINDEX_END'};
  $grid_msg = "I am $grid_mode job $jobid task $taskid of $maxtask \@ $ENV{HOSTNAME}\n";
}

# Debugging
warn $grid_msg if defined $grid_msg && $CONFIG{'GRID_ARRAY'} && !$opts{'q'};

# If -q option was specified then query file has already been processed in PPH::Config
my $input_file;
unless ($opts{'q'}) {
  for (my $i=0; $i<@ARGV; $i++) {
    my $abspath = abs_path($ARGV[$i]);
    die "ERROR: Input file does not exist: $ARGV[$i]\n"
      unless defined $abspath && length $abspath && -e $abspath;
    $ARGV[$i] = $abspath;
  }
  # Remember the name of the first input file (if any) before it gets shifted out by <> below
  $input_file = $ARGV[0];
}

# Set current stripe index and total number of stripes
my ($stripe_no, $stripes);
if  ($opts{'r'}) {
  if      ($opts{'q'}) {
    warn "WARNING: Input stripe specification via '-r' option ignored in query mode\n";
  } elsif ($array_mode) {
    warn "WARNING: Explicit input stripe specification via '-r' option ignored in $grid_mode array mode\n";
  } else {
    ($stripe_no, $stripes) = $opts{'r'} =~ m|^(\d+)/(\d+)$|;
    die "ERROR: Illegal option's argument: -r $opts{'r'}\n" unless $stripe_no && $stripes;
  }
}
if ($array_mode) {
  ($stripe_no, $stripes) = ($taskid, $maxtask);
}

# Calculate range of input lines if striping option was specified
my $input_line_count  = 0;
my $input_query_count = 0;
my ($range_start, $range_end);
if ($stripes) {
  # Count how many input lines we have
  for (my $i=0; $i<@ARGV; $i++) {
    open(FIN, $ARGV[$i]) or die "ERROR: Can't open input file: $ARGV[$i]\n";
    while (<FIN>) {
      $input_line_count++;
      # Skip comments and empty lines
      next if /^#/ || /^\s*$/;
      $input_query_count++;
    }
    close(FIN);
  }
  # Calculate stripe range
  ($range_start, $range_end) = stripe_range($stripe_no, $stripes, $input_query_count);
  # Validate stripe range
  if ($range_start == 0 && $range_end == 0) {
    # In array mode, if $maxtask was set to a value higher than the total number of queries
    # in input then we kill the rest of array as soon as the current task gets beyond the
    # input range, to avoid wasting cluster slots on keeping validating out-of-range tasks.
    # Note: range_start/range_end are both zero if and only if stripe_no > input_query_count
    if ($array_mode) {
      if ($maxtask) {
        if      ($grid_mode eq 'SGE') {
          warn "WARNING: SGE task $jobid.$taskid exceeds input range, rest of array omitted\n";
          if ($stripe_no < $maxtask) {
            my $qdel = "sgesh_anko qdel $jobid -t " . ($stripe_no + 1) . "-$maxtask";
            my $qrc  = `$qdel 2>&1`;
            warn "WARNING: qdel error: $qrc\n" unless
              $qrc =~ /^\w+ has deleted job-array tasks \d+-\d+:\d+ of job \d+/ ||
              $qrc =~ /^\w+ has registered the job-array task \d+\.\d+ for deletion/;
          }
        } elsif ($grid_mode eq 'LSF') {
          warn "WARNING: LSF task $jobid.$taskid exceeds input range, rest of array omitted\n";
          if ($stripe_no < $maxtask) {
            my $qdel = "bkill -b '${jobid}[" . ($stripe_no + 1) . "-${maxtask}]'";
            my $qrc  = `$qdel 2>&1`;
            warn "WARNING: bkill error: $qrc\n" unless
              $qrc =~ /^The requested operation is in progress/ ||
              $qrc =~ /^Job has already finished/;
          }
        } else {
          warn "WARNING: Task $jobid.$taskid exceeds input range, terminated\n";
        }
        exit 0;
      }
    }
    die "ERROR: Input out of range (stripe $stripe_no / $stripes with $input_query_count queries in input)\n";
  }
}

my $PPH = PPH->new($CONFIG{'SEQFILE'});

# Create a temporary directory inside TMP_OUT and chdir into it.
# Also set TMP_OUT to point to this directory. Save old CWD in
# case we later need to know where it pointed originally.
$CONFIG{'TMP_OUT'} = tempdir("$CONFIG{TMP_OUT}/PPH_XXXXXXXX");
$CONFIG{'CWD'} = cwd();
chdir $CONFIG{'TMP_OUT'} or die "WARNING: Failed to chdir into temporary directory: $CONFIG{TMP_OUT}\n";
### |   Scratch directory is set to    : $CONFIG{SCRATCH}
### |   Temporary directory created as : $CONFIG{TMP_OUT}

# Print header (formatted) unless running in array mode
print('#', join("\t", $PPH->return_header(-formatted=>1)), "\n") unless $array_mode || $opts{'r'} && $stripe_no != 1;

#----------------------------------------
# Read input
#----------------------------------------
my $input_line_no  = 0;
my $input_query_no = 0;
my $range_line_no  = 0;
my $query_no = 1;
my $input_no = 0;
while ($_ = $CONFIG{'QUERYSTRING'} || <>) {

  $input_line_no++;   # unfiltered input line number, used in error/warning messages
  $input_no++;

  # Skip comments and empty lines
  next if /^#/ || /^\s*$/;

  $input_query_no++;  # input query line number

  if ($stripes) {
    next if $input_query_no < $range_start;
    last if $input_query_no > $range_end;
  }

  $range_line_no++;   # total query input lines to process (i.e., within range but without comments)

  chomp;

  # Helps debugging of parallel tasks executed on different cluster nodes
  warn "INPUT: $_\n" if $array_mode;

  # Trim optional user comments
  my $comments = $1 if s/\s *#\s*(.*)$//;

  # In input, we can have either a single dbSNP rsID column, or a 4-column mandatory
  # AA residue substitution specification (accession, position, aa1, aa1) optionally
  # followed by up to 10 columns of nucleotide parameters (txname, cdspos, frame,
  # nt1, nt2, flanks, transv, cpg, jxdon, jxacc).
  # Note: we split first 4 mandatory columns on whitespace but then split the rest on
  #       tabs to preserve any empty optional fields (e.g., jxdon and jxacc).
  my @params;
  my $paracount = 0;
  # Mandatory parameters
  my @mandatory = split ' ', $_, 5;
  $paracount    = scalar @mandatory;
  if      ($paracount == 4 || $paracount == 1 && $mandatory[0] =~ /^rs\d+$/) {
    @params = @mandatory;
  } elsif ($paracount < 4) {
    warn "ERROR: Insufficient number of parameters ($paracount) in input at line $input_no\n";
    last if $CONFIG{'QUERYSTRING'};
    next;
  # Optional parameters present
  } elsif ($paracount > 4) {
    my @optional  = split /\t/, $mandatory[-1], 11;
    my $present   = 0;
    map { s/^\s+//; s/\s+$//; $present++ if length } @optional;
    if ($present > 10) {
      warn "ERROR: Too many parameters ($present) in input at line $input_no\n";
      last if $CONFIG{'QUERYSTRING'};
      next;
    }
    @params = (@mandatory[0..($#mandatory-1)], @optional);
  }

  #------------------------------------------
  # Creating and verifying new SNP object
  #------------------------------------------
  my $SNP;
  eval { $SNP = SNP->new(\@params, $comments) };
  if ($@) {
    warn "ERROR: Verification failed for input line $input_no: $@";
    last if $CONFIG{'QUERYSTRING'};
    next;
  }

  # Failure of prepararation if non-empty msg returned
  my $msg = $PPH->prepareInput($SNP);
  if ($msg) {
      warn "ERROR: $msg\n";
      last if $CONFIG{'QUERYSTRING'};
      next;
  }
  #------------------------------------------

  #------------------------------------------
  # Running PolyPhen-2
  #------------------------------------------
  $PPH->analyze($SNP);
  #------------------------------------------

  #------------------------------------------
  # Evaluating data
  #------------------------------------------
  $PPH->interpret_results($SNP);
  #------------------------------------------

  # Wrap current SNP object into PPH object
  push @{ $PPH->{Variant} }, $SNP;

  #-----------------------------------------------
  # Printing tab-delimited summary table to STDOUT
  #-----------------------------------------------
  print join("\t", $PPH->format_result($SNP)), "\n";
  #-----------------------------------------------

  last if $CONFIG{'QUERYSTRING'};
}

#-----------------------------------------------------
# Save PPH object to a file in XML format (optional,
# currently only used by the web interface).
#-----------------------------------------------------
if ($CONFIG{'XMLOUTPUT'}) {
  # Add sources versioning information if available
  $PPH->{Source}{Predictor} = "PolyPhen-2 v$VERSION" if defined $VERSION;
  if (exists $CONFIG{'UNIPROT'}) {
    my $version_file = dirname($CONFIG{'UNIPROT'}) . '/VERSION';
    if (-s $version_file) {
      my $desc;
      if (open my $F, '<', $version_file) {
        chomp($desc = <$F>);
        close $F;
      }
      $PPH->{Source}{Sequence} = $desc if defined $desc && length $desc;
    }
  }
  if (exists $CONFIG{'PDB2FASTA'}) {
    my $version_file = $CONFIG{'PDB2FASTA'} . '/VERSION';
    if (-s $version_file) {
      my $desc;
      if (open my $F, '<', $version_file) {
        chomp($desc = <$F>);
        close $F;
      }
      $PPH->{Source}{Structure} = $desc if defined $desc && length $desc;
    }
  }
  if (exists $CONFIG{'MULTIZPATH'}) {
    my $version_file = $CONFIG{'MULTIZPATH'} . '/precomputed/VERSION';
    if (-s $version_file) {
      my $desc;
      if (open my $F, '<', $version_file) {
        chomp($desc = <$F>);
        close $F;
      }
      $PPH->{Source}{Mz} = $desc if defined $desc && length $desc;
    }
  }
  my $xmlfile;
  if ($array_mode) {
    $xmlfile = "$jobid.$taskid";
  } else {
    if      ($opts{'q'}) {
      if ($grid_mode) {   # hack for the GGI single query mode
        $xmlfile = "$jobid";
      } else {
        $xmlfile = fileparse($opts{'q'});
        $xmlfile =~ s/\.[^.]*$//;
      }
    } elsif (defined $input_file) {
      $xmlfile = fileparse($input_file);
      $xmlfile =~ s/\.[^.]*$//;
      $xmlfile =~ s/\.pph$//;
    } else {
      $xmlfile = 'stdin.' . $$;
    }
  }
  # Note: use saved CWD for output target directory path since
  # the script may be run chdir()'ed into a temporary directory
  $xmlfile = $CONFIG{'CWD'} . '/' . $xmlfile . '.pph.xml';
  open my $F, '>', $xmlfile or die "ERROR: Can't create output file: $xmlfile\n";
  $PPH->pph2xml($F);
  close $F;
}
#-----------------------------------------------------

warn "TOTAL: $range_line_no queries processed\n" if $array_mode;

#----------------------------------------------------------------------------------------

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

END {
  # Newer rmtree() versions refuse to remove directory which is
  # still current so we need to first go back to where we belong
  chdir $CONFIG{'CWD'} if exists $CONFIG{'CWD'};
  if (!$DEBUG && exists $CONFIG{'CWD'} && -d $CONFIG{'TMP_OUT'}) {
    rmtree($CONFIG{'TMP_OUT'})
      or warn "Failed to remove temporary directory: $CONFIG{TMP_OUT}\n";
  }
}
