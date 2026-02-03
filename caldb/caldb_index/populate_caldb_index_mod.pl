#!/usr/bin/perl
#
# Given an CALDB index file and a root level directory, traverse the directory structure to include
# all subdirectories for a complete list of all files, then update those to the specified CALDB
# index file using FTOOLS command 'udcif'
#
# Uses Perl File::Spec library for tracking path information
#
# 21-Aug-2024: Added flag to write commands to a file without running them during the
#              execution of the script. (James Runge)

use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my $VERSION = '0.1 (21 August 2024)';

# Requires FTOOLS and CALDB setup: exit is not configured for FTOOLS

unless (exists $ENV{HEADAS} && exists $ENV{CALDB}) {
    print "Requires FTOOLS to be setup and CALDB to be defined\n";
    exit;
}

my $help = 0;
my $verbose = 0;
my $root_dir = '';
my $caldb_index = '';
my $exclude_dirs = '';
my $outfile = 0;
my $filename = '';
my $fh;
my $result = GetOptions( "rootdir=s" => \$root_dir,
			 "index=s"   => \$caldb_index,
			 "exclude=s" => \$exclude_dirs,
			 "outfile" => \$outfile, 
                         "verbose"  => \$verbose,
                         "help"     => \$help );
if (! $result)  { use_message(); exit; }
if ($help)      { use_message(); exit; }
if ($#ARGV > 0) { 
    print "Command line argument not recognized\n"; 
    use_message(); 
    exit; 
}
# Error conditions
unless (-d $root_dir) {
    print "Unable to find root level directory $root_dir\n";
    exit;
}

if ($root_dir !~ /\/$/) { $root_dir .= '/'; }
unless (-f "${root_dir}${caldb_index}") {
    print "Unable to find CALDB index file $caldb_index\n";
    exit;
}

# If outfile flag, write commands to a file
if ($outfile) {
    $filename = 'caldb_commands.csh';
    open($fh, '>', $filename) or die "Could not open file '$filename' $!";
}

# Exclude the directory containing the CALBD index file from consideration
# then also exclude any user-requested directories to exclude

my ($volume, $directory, $file) = File::Spec->splitpath($caldb_index);
my @excludes = ( $directory );
if ($exclude_dirs) {
    my @dirs = split /\,/, $exclude_dirs;
    # Make sure each excluded dir as a trailing '/'
    foreach my $dir (@dirs) {
	unless ($dir =~ /\/$/) { $dir .= '/'; }
        push @excludes, $dir;
    }
}

# Command to execute on each file, maintaining relative path to cif (CALDB index file)

# Find and process each file found under $root_dir
if ($verbose) { print "Finding all items in $root_dir\n"; }
my @dirtree = `find $root_dir`;
if ($? != 0) { print "Failure when using find command on $root_dir\n"; exit; }

# Prune the directory tree to make an array that only lists files and includes
# no files from any of the excluded directories in the directory tree
if ($verbose) { print "Pruning out directories and files in excluded directories\n"; }
my @files = filter_directories(\@dirtree,\@excludes);

if (scalar @files == 0) {
    print "Failed to find any files\n";
    exit;
}
if ($verbose) {
    printf "Found %d files in %s\n", scalar @files, $root_dir;
}

# Select all the files in the pruned directory tree and create index entries for each
my $prior_directory = '';
foreach my $filename (@files) {
    # Split out the full path and the file name
    ($volume, $directory, $file) = File::Spec->splitpath($filename);
    if ($directory ne $prior_directory) {
	if ($verbose) { print "Moving to $directory\n"; }
	if ($outfile) { print $fh "cd $directory\n"; }
	chdir $directory;
	$prior_directory = $directory;
    }
    # Convert the $caldb_index setting into a relative file name
    my $relative_index = File::Spec->abs2rel($root_dir . $caldb_index, $directory);
    my $command = 'udcif infile="' . $file . '" cif="' . $relative_index . '" quality=0';
    my $len = length($file);
    if ($verbose) {
	print "$command\n";
    }
    if ($outfile) {
	print $fh "echo '$file'\n";
	print $fh "$command\n";
    }
    else {
	my $status = system($command);
	if ($status != 0) {
	    # Error processing
	    print "Error running $command\n";
	    exit;
	}
    }
}


sub filter_directories {
    my ($tree_ref, $exclude_ref) = @_;

    my @files_found = ();
    my @entries = @{$tree_ref};
    my @exclude_dirs = @{$exclude_ref};

    foreach my $item (sort @entries) {
	chomp $item; # Trailing \n from using of backticked find command
	my $excluded = 0;
	unless (-f $item) { next; }
	foreach my $exclude_dir (@exclude_dirs) {
	    if ($item=~ /^${root_dir}${exclude_dir}/) {
		$excluded = 1;
		last;
	    }
	}
	push @files_found, $item unless $excluded;
    }

    return @files_found;
}


# Standard Usage message

sub use_message {
    print <<"EndOfMessage"
VERSION: $VERSION

populate_caldb_index.pl --rootdir full_path --index index_file [--exclude ex_dirs] [--outfile] [--verbose] [--help]

where  full_path    is the root level directory to start the search for all files
       index_file   is the location and name of the CALDB index file 
                    relative to full_path
       ex_dirs      is an optional comma-separated list of directories under full_path to not
                    include in the CALDB index
       outfile      is an optional flag to write commands to a file to be run later without
                    executing them at runtime
       verbose      is an optional flag to run the code in a more talkative mode
       help         prints this help message

Update an existing/empty CALDB index file "index_file" with all the files found below the directory
 full_path, and adding each such file to "index_file" using udcif. Exclude all files in the 
 same path as the index file. If requested, also exclude any files in any of a list of excluded
 directories.

NOTE:
Requires CALDB and HEADAS environments to be set
Requires CALDB to be set to an editable version of CALDB
The command udcif will produce output, so even with --verbose not enabled, there may be considerable
 output if there are a large number of files to be indexed.


EXAMPLE

% populate_caldb_index.pl --rootdir /processing/calet/caldb/data/calet/cgbm --index index/caldb.indx_tmp --exclude create_index --verbose

This updates caldb.indx_tmp to contain entries for every single file found below 
/processing/calet/caldb/data/calet/cgbm EXCEPT any files in index/ (so CALBD goes not get 
entangled in itself) nor any files in the directories create_index/ (so, in this example, the 
script populate_caldb_index.pl and other items unrelated to CALDB could be in create_index/
and not be included.

This is intended to build an entire new CALDB mission entry in a single pass with a single command.
This is especially useful for calibration databases that have a very large number of entries.

EndOfMessage
}
