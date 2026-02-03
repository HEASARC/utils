#!/usr/bin/perl

# NAME: checksum_transfer.pl
#
# AUTHOR: David Riethmiller, NASA GSFC
#
# PURPOSE: Verify that a data set has not been truncated accidentally during transfer, by 
#          comparing starting and finishing checksums.
#
# CALLING SEQUENCE:  checksum_transfer.pl [TOPDIR FILE] [LOCATION] [ORIGIN FILE] [DESTINATION FILE]
#
# INPUT:  TOPDIR FILE: ascii file listing full paths to all top directories to be checked.
#         LOCATION: location at which checksum is being done; "O" for origin, "D" for destination
#         ORIGIN FILE: name of file containing checksums from origin location (Origin mode: will be 
#           written, Destination mode: will be read)
#         DESTINATION FILE: name of file containing checksums from destination location (only used 
#           in Destination mode)
#         *** these parameters will be prompted for if not provided.
#
# OUTPUT:  an ascii text file listing recursively every file under the input top directory list and 
#          its checksum values.  If doing a destination check, this file is compared against the 
#          specified origins checksum file.
#
# Revised 2019-11-19 to include sorting on filenames to erase file listing order differences 
# as a source of origin/destination mis-matches. Added strict and warnings for Perl stability.
# J. Allen NASA/GSFC/HEASARC
# Revised 2019-11-25 to tolerate trailing spaces and/or blank lines in directory listings

use strict;
use warnings;
use File::Copy qw(move);
use File::Find;

my $dirfile = "";
my $org_or_dest = "";
my $origin_file = "";
my $destination_file = "";

# Get the list of topdirs included in the data transfer.
if ( defined $ARGV[0] ){
    $dirfile=$ARGV[0];
} else {
    print "\nEnter name of ascii file listing full path name to every top directory to be checked (i.e., dirlist.txt):\n%> ";
    $dirfile=<STDIN>;
}
chomp $dirfile;

print "\ndirfile = $dirfile\n";

# Ask the user whether this is an origin or destination check.
my $valid_OD = 0;
while ($valid_OD == 0){

    if ( defined $ARGV[1] ){
	$org_or_dest = $ARGV[1];
    } else {
	print "\nIs this an Origin (letter O) or Destination (letter D) check?\n%> ";
	$org_or_dest = <STDIN>;
    }
    chomp $org_or_dest;
    $org_or_dest = uc $org_or_dest;

    if ("$org_or_dest" eq "O") {
	$valid_OD = 1;
	print "\nDoing ORIGIN check.\n";
    } elsif  ("$org_or_dest" eq "D") {
	$valid_OD = 1;
	print "\nDoing DESTINATION check.\n";
    } else {
	print "That is not a valid location.  Please enter letters O or D.\n";
    }
}

# Get the name of the origin/destination files
if ("$org_or_dest" eq "O") {

    if ( defined $ARGV[2] ){
	$origin_file = $ARGV[2];
    } else {
	print "\nEnter name of origin file to be written (i.e. checksums_origin.txt):\n%> ";
	$origin_file=<STDIN>;
    }
    chomp $origin_file;

    if ($origin_file eq ""){ # Make sure the user actually entered something
	print "   ERROR: must enter origin file name.\n";
	exit;
    } else {
	print "\nOrigin file: $origin_file \n";
    }

} elsif ("$org_or_dest" eq "D") {

    if ( defined $ARGV[2] ){
        $origin_file = $ARGV[2];
    } else {
	print "\nEnter name of origin file to be checked against (i.e. checksums_origin.txt):\n%> ";
	$origin_file=<STDIN>;
    }
    chomp $origin_file;

    if (-e $origin_file) {
	print "\nLocated file $origin_file successfully.\n";
    } else {
	print "   ERROR: could not find file $origin_file.\n";
	exit
    }

    if ( defined $ARGV[3] ){
        $destination_file = $ARGV[3];
    } else {
	print "\nEnter name of destination file to be written (i.e. checksums_destination.txt):\n%> ";
	$destination_file=<STDIN>;
    }
    chomp $destination_file;

    if ($destination_file eq ""){ # Make sure the user actually entered something
        print "   ERROR: must enter destination file name.\n";
        exit;
    } else {
	print "\nDestination file: $destination_file\n";
    }

} else {
    print "This condition should never happen.\n";
    exit;
}

# Read the topdir file and populate the topdir list
my @topdir_full;  # The full path to the top directory listed, i.e. /highest/high/med/low/lowest
my @topdir_base;  # Only the last field of the top directory path, i.e. /lowest
open (my $FILE1, $dirfile);
while (my $row = <$FILE1>){
    chomp $row;
    trim($row);
    if (length($row) > 0){
        $row =~ s/\s+$//; # Trim out any trailing blank space
	push @topdir_full,$row;
    } else { next; } # Skip over blank lines
    if (! -d $row){
	print "ERROR: could not find directory $row, exiting.\n";
	exit;
    }
    # Isolate the topdir number and save to array
    my $obs_only = `basename $row`;
    push @topdir_base, $obs_only
}
close $FILE1;


my $tmp_out = "tmp_checksums.txt";
chomp $tmp_out;

# Open a temporary output file for writing.  We'll rename it later.
open( my $FILE2, '>', $tmp_out);

# Cycle through the topdir directories
my $dir_index = 0;
foreach my $dir (@topdir_full) {
    my @allfiles;
    my $dirnum = $topdir_base[$dir_index];
    $dir_index++;
    

    # Get all the files under each topdir 
    find sub {
	return if -d;
	push @allfiles, $File::Find::name;
    }, $dir;

    # Compute a checksum for every file and write to log
    foreach my $file (@allfiles){
	my $cksum=`cksum $file | awk -F ' ' '{print \$1}'`;
	chomp $cksum;

	# Cut everything to the left of the lowest topdir from the file name
	my $len_dirnum = length $dirnum;
	my $len_dir = length $dir;
	my $NN = $len_dir - $len_dirnum;
	my $to_remove = substr($file, 0, $NN) . "/";
	$file =~ s/$to_remove//;
	chomp $file;

	#printf "%-15s  %s\n", $cksum, $file;
	printf $FILE2 "%-15s  %s\n", $cksum, $file;
    }
}

close $FILE2;


print "\n";

if ("$org_or_dest" eq "D")  { # If this is the destination, write the desination file and compare against origin file
    move $tmp_out, $destination_file;
    print "Comparing against origin checksums...\n\n";
#    my $diff=`diff $origin_file $destination_file`;
# 2019-11-19 Sort both files to remove file listing order differences unrelated to checksums
    print "Making sorted temporary version of $origin_file\n";
    if (-e "tmp_origin_sorted.txt") { unlink("tmp_origin_sorted.txt"); }
    if (-e "tmp_destination_sorted.txt") { unlink("tmp_destination_sorted.txt"); }
    my $status = system("sort -k 2 $origin_file > tmp_origin_sorted.txt");
    if ($status != 0) 
        { print "WARNING: sort failed to map $origin_file to tmp_origin_sorted.txt\n"; }
    print "Making sorted temporary version of $destination_file\n";
    $status = system("sort -k 2 $destination_file > tmp_destination_sorted.txt");
    if ($status != 0) 
        { print "WARNING: sort failed to map $destination_file to tmp_destination_sorted.txt\n"; }
    my $diff=`diff tmp_origin_sorted.txt tmp_destination_sorted.txt`;
    if ( $diff eq ""){
	print "Files $origin_file and $destination_file are identical; all checksums verified.\n";
        unlink("tmp_destination_sorted.txt");
        unlink("tmp_origin_sorted.txt");
    } else {
	print $diff;
    }

} elsif ("$org_or_dest" eq "O") { # If this is the origin, write the origin file
    print "Creating origin checksums file $origin_file.\n";
    move $tmp_out, $origin_file;
}

print "\n";


# End of main program
exit;



# Subroutine to trim white space.
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
	
