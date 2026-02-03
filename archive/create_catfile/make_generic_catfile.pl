#!/usr/bin/perl
#
# Read a generic template to read an archive or subsection of an archive and create matching
# FITS catalog (.cat) files
# 
# Author: Jesse Allen (NASA/GSFC HEASARC)
#
# VERSION 0.801
#
# V 0.1  2019-11-02  
# V 0.2  2019-12-08
#  Added "range" setting to permit limits on search for variables
# V 0.3  2020-02-10
#  Removed automatic auxil/ path for catfile, must instead specify where cat goes
#  in template file (to accomodate HaloSat, which has no auxil directory)
# V 0.4  2020-03-02
#  Replaced use of maxitool for mjd2date with computing approximate unixtime
#   (seconds since 1970-01-01T00:00:00, ignoring leap seconds)
#   and using strftime POSIX function to generate year, month, day from MJD
#   Add CATDEST (path within the archive to the final archive .cat file) to the
#    catalog line entry
# v 0.5  2021-09-30
#   Added --range (variable) --min (min_value) --max (max_value) to command line
#   Allows ONE settings only to be limited from the commmand line
# v 0.6  2023-06-29
#   Interpret wildcards in filelist subdirectories in full in order to get files 
#   more than a single directory level below the beginning of the observation directory
# v 0.7  2023-07-11
#   Get keywords from first FITS file in list of files (previously got from first file
#   without checking if FITS)
# v 0.8  2025-08-08
#   Improve some error trapping
#    Report on missing path
#    Clean path if trailing whitespace
#    Report failure to find a keyword
#   Use ftkeypar instead of fkeypar (should both work the same, but may not...)
#   Skip over existing .cat or .cat.gz file
#   Support use of variables in categories
#     Trap error condition if variable in category is not known
# v 0.801  2025-08-15
#   Check if keyword exists for TUNIT cleanup step
#
# TO DO
# Check on OBSID vs OBS_ID keyword
#
# Code currently has "baked in" format for CAT file, not fully attuned to cdfile 
# settings: look at ASCII template print statement and adjust
#
# WARNINGS/CONCERNS
# .cat files can be compressed, but the field settings in the .cat will ALWAYS show
# only uncompressed settings

use strict;
use warnings;
use POSIX;
use Getopt::Long;  # Get Command line options for template file
use Astro::FITS::CFITSIO qw(:longnames :constants);

my $CODE    = 'make_generic_catfile.pl'; # This program
my $VERSION = '0.8 (8 Aug 2025)';        # Version and date of this program

my $mjdref = 40857; # MJD on 1970-01-01 (unixtime 0)
my $template;
my $cdfile;
my $primary;
my $extension;
my $help = 0;
my $verbose;
my $compress_cat = 0;
my $archive_cat = 0;
my $range_setting = undef;
my $range_min = undef;
my $range_max = undef;

GetOptions( "template=s"  => \$template,
            "cdfile=s"    => \$cdfile,
            "primary=s"   => \$primary,
            "extension=s" => \$extension,
            "archive"     => \$archive_cat,
            "range=s"     => \$range_setting,
            "min=i",      => \$range_min,
            "max=i"       => \$range_max,
            "compress"    => \$compress_cat,
            "help"        => \$help,
            "verbose"     => \$verbose )
    or use_message(); 

if ($help) { use_message(); exit; }
if (!-e "$template")  { die "Could not find $template"; }
if ($cdfile eq '') { 
    $cdfile = 'cdfile_basic.txt';
    print "No cd file (columns and settings) provided, using default CATFILE settings\n";
}
unless (-e "$cdfile") 
    { die "Could not find cdfile $cdfile"; }
unless (-e "$primary") 
    { die "Could not find primary header keywords template file $primary"; }
unless (-e "$extension") 
    { die "Could not find 1st extension header keywords template file $extension"; }
if ($archive_cat) { 
    print "Once catalog files are created, they will be placed in the archive\n"; 
} else {
    print "Once catalog files are created, they will remain in the local directory\n"; 
}
if (defined $range_setting) {
    unless ((defined $range_min) && (defined $range_max)) {
        die "Must specify --min and --max when used with --range";
    }
}
print "$CODE $VERSION\n";

# Set hash tables of FITS keywords to be used in the primary array and first
# extension of the CAT file 
my ($fixed_ref, $filename_ref, $updated_ref) = read_keyword_template($primary);
my %primary_fixed_keywords = %{$fixed_ref};
my %primary_filename_keywords = %{$filename_ref};
my %primary_updated_keywords = %{$updated_ref};

($fixed_ref, $filename_ref, $updated_ref) = read_keyword_template($extension);
my %firstext_fixed_keywords = %{$fixed_ref};
my %firstext_filename_keywords = %{$filename_ref};
my %firstext_updated_keywords = %{$updated_ref};

# Open and read the template file decribing which files in the archive should be 
# cataloged, their characterstics, and the pattern of filenames and paths to them.
# These include "variables" that may be in filenames and paths, such as year and month,
# MJD date, and observation ID.
open (my $tmp_fh, '<', $template);
my $counter = 0;
# Patterns use Unix-type wildcards, not regex or Perl regex 
# since file searching uses "find" to discover files in the archive.
# %pattern is initialized with some of the most common such settings, but additional
# patterns may be specified in the template, or existing ones may be overridden with 
# differing patterns 
my %pattern = (
    'YYYY',  '[12][0-9][0-9][0-9]',  # Years
    'YY',    '[0-9][0-9]',           # Years without leading century (discouraged)
    'MM',    '[01][0-9]',            # Months with leading 0
    'DD',    '[0-3][0-9]',           # Days of the month with leading 0
    'DDD',   '[0-3][0-9][0-9]',      # Julian day of the year (1-365 (366 if leap year)
    'OBSID', '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]',
                                     # Default OBSID format is 10-digit number
    'MJD1000', '[0-9][0-9]000',              # MJD 1000 is int(MJD/1000) * 1000
    'MJD',     '[0-9][0-9][0-9][0-9][0-9]'   # MJD is always a five digit number
    );
my $path;
my $subdir;
my $catname;
my $catdest;
my $filename;
my @filelist_variables;
my $category_counter;
my %filetype;
my %var_min;
my %var_max;

if (defined $range_setting) {
    $var_min{$range_setting} = $range_min; 
    $var_max{$range_setting} = $range_max;
}
while (<$tmp_fh>) {
    if (/\#/) { next; } # Skip comment lines
    chomp;
    # Read in settings
    if (/^variable (ROOT)?PATH (.+)$/) { $path = $2; next; }
    if (/^variable SUBDIR (.+)$/)      { $subdir = $1; next; }
    if (/^variable CAT(DEST|DIR) (.+)$/) { $catdest = $2; $pattern{'CATDEST'} = $2; next; }
    if (/^variable CATNAME (.+)$/)     { $catname = $1; $pattern{'CATNAME'} = $1; next; }
    # Allows user-specific variable or overriding defaults (above)
    if (/^variable (\w+) (.+)$/) {
        $pattern{$1} = $2;
        next;
    }
    # Allows range limits to be set on variables, overriden by command line
    if (/^range (\w+) (.+) (.+)$/) {
        if ((defined $range_setting) && ($range_setting eq $1)) {
            if ($verbose) { 
                print "RANGE $1 overridden by command line setting\n";
	    }
	} elsif ($2 lt $3) {
            $var_min{$1} = $2;
            $var_max{$1} = $3;
	} else {
            $var_min{$1} = $3;
            $var_max{$1} = $2;
	}
        next;
    }
    # Read in list (or perhaps several lists) of files to be listed in CATFILE
    if (/^filelist/) {
        my @elements = split;
        $filelist_variables[$counter++] = $elements[1];
        next;
    }
    # Read in list (or perhaps several lists) of files to be listed in CATFILE
    if (/^category/) {
        my @elements = split(/\s+/,$_,6);
        my $key = $elements[1];
        my $descript = "'" . $elements[5] . "'";
        $filetype{$key} = [$elements[2], $elements[3], $elements[4], $descript];
        next;
    }
    # Error checking. Lines should start with # (comment), range, category, variable,
    # or filelist or blank line 
    if (!/^(range|category|variable|filelist|\#|\s*$)/) {
	print "WARNING! Line does not have expected contents: $_\n";
    }
}
close $tmp_fh;
if ($counter == 0) { die "No file list specifications found in template $template"; }

$path =~ s/\s*$//;   # Remove any trailing space
$subdir =~ s/\s*$//; # Remove any trailing space
# Add a trailing '/' to path and subdir if necessary
# Exit is path to files does not exist
if ($path) {
    if (substr($path,length($path)-1,1) !~ /\//) { $path .= '/'; }
    if (!-d $path) {
	print "Could not locate $path, exiting\n";
	exit;
    }
}
if ($subdir) {
    if (substr($subdir,length($subdir)-1,1) !~ /\//) { $subdir .= '/'; }
    # Do NOT test for existence of subdirectory since it may have an unassigned variable
    #if (!-d "${path}${subdir}") {
    #	print "Could not locate ${path}${subdir}, exiting\n";
    #	exit;
    #}
}

# Parse the variables, if any, in the subdirectory path into pattern matches
# for glob
my $path_search;
my @vars_in_path;
my @patterns_in_path;
if ($subdir =~ /:\w+:/) {
    ($path_search,my $var_ref,my $pattern_ref) = replace_varname_with_pattern($path.$subdir, 
                                                                              \%pattern);
    @vars_in_path = @{$var_ref};
    @patterns_in_path = @{$pattern_ref};
} else {
    print "No variables in subdir or no subdir provided\n";
    $path_search = $path . $subdir;
}
if ($path_search =~ /\s+$/) {
    print "Found trailing space in path search, removing...\n";
    $path_search =~ s/\s*$//;
}

# Create sorted array of file type suffixes, with longest suffixes first (to avoid confusion
# between, for example, ".att.hk" and ".hk".)
# Check for variables embeded in filetype categories and replace them with patterns

my %new_filetype;
my @unsorted;
my $i = 0;
for my $key (keys %filetype) {
    my ($category, $var_ref, $pattern_ref) =  replace_varname_with_pattern($key, \%pattern);
    $new_filetype{$category} = $filetype{$key};
    $unsorted[$i++] = $category;
}
# Update %filetype with the pattern revisions
%filetype = %new_filetype;
my @sorted = sort { length($b) <=> length($a) } @unsorted;

# Parse the variables, if any, in the file lists
my @file_search_list;
my @vars_in_filelist;
my @patterns_in_filelist;
for (my $i = 0; $i < $counter; $i++) {
    if ($filelist_variables[$i] =~ /:\w+:/) {
        my ($file_search,$var_ref,$pattern_ref) = 
            replace_varname_with_pattern($filelist_variables[$i],\%pattern);
        $file_search_list[$i] = $file_search;
#        $vars_in_filelist[$i] = @{$var_ref};
#        $patterns_in_filelist[$i] = @{$pattern_ref};
    } else {
        print "No variables in file list provided\n";
        $file_search_list[$i] = $filelist_variables[$i];
    }
}

# Search for all the variations of path + subpath
# Limit how deep find looks (or it might try to load the entire archive file listing!)
my @path_array = `find $path_search -maxdepth 0`;
for my $all_path (@path_array) {
    chomp($all_path); 
    if ($all_path !~ /\/$/) { $all_path .= '/'; }
    my %current_value = find_value_with_pattern($path.$subdir, 
                                                $all_path, 
                                                \@vars_in_path, \@patterns_in_path);

# Translate CATNAME to a complete filename. Note this means CATNAME must 
# have a variable name from the path

    my $file_pattern = $pattern{'CATNAME'}; 
    my $yr; my $mn; my $dy; my $mjd_date; my $date_obs; my $obsid_date;
    my $out_of_range = 0;
    if ($file_pattern =~ /\:\w+\:/) {
	foreach my $var_to_check (@vars_in_path) {
            my $current_setting;
            if ($file_pattern =~ /\:${var_to_check}\:/) {
                $current_setting = $current_value{$var_to_check};
                if ($file_pattern =~ /:YYYY:/) { $yr = $current_setting; }
                if ($file_pattern =~ /:YY:/) { 
                    if ($current_setting < 69) { $yr = 2000 + $current_setting; }
                    else                       { $yr = 1900 + $current_setting; }
	        }
                if ($file_pattern =~ /:MM:/)   { $mn = $current_setting; }
                if ($file_pattern =~ /:DD:/)   { $dy = $current_setting; }
                if ($file_pattern =~ /:MJD:/)  { $mjd_date = $current_setting; }
		$file_pattern =~ s/\:${var_to_check}\:/${current_setting}/;
	    } elsif (($file_pattern =~ /:YY:/) && ($var_to_check eq 'YYYY')) {
                $current_setting = substr($current_value{$var_to_check},2,2);
                $yr = $current_value{$var_to_check};
		$file_pattern =~ s/\:${var_to_check}\:/${current_setting}/;
	    }
            # Set flag to indicate to skip this file if any variable is out of user-specified range
            if (exists ($var_min{$var_to_check})) {
                if (($current_setting lt $var_min{$var_to_check}) ||
                    ($current_setting gt $var_max{$var_to_check})) {
                    $out_of_range = 1;
	        }
	    }
	}
        if ($out_of_range) { next; } # Skip to next file if any variable is out of range
	if ($mjd_date) { 
            my $unixtime = ($mjd_date - $mjdref) * 86400;
            $yr = strftime("%Y", $unixtime);
            $mn = strftime("%m", $unixtime);
            $dy = strftime("%d", $unixtime);
	    $obsid_date = $mjd_date;
	} elsif ($yr) {
           $obsid_date = sprintf("%04d%02d%02d", $yr, $mn, $dy);
	}
        if ($yr) { $date_obs = sprintf("%04d-%02d-%02d", $yr, $mn, $dy); }
        # Now fail of there are any variables in the filename that could not be derived from path
        if ($file_pattern =~ /\:(\w+)\:/) { die "Could not match $1 from variables in path"; }
    }
    my $catfile = $file_pattern;
    my @full_catlist;
    # Now create a list of all files to be placed in the CATFILE
    for (my $i = 0; $i < $counter; $i++) {
        my @dirs = split(/\//, $all_path . $file_search_list[$i]);
        pop @dirs;
        my $testpath = join('/',@dirs);
        # If testpath has wildcards, interpret them
        my @checkdirs = glob "$testpath";
        if (scalar @checkdirs == 0) {
	    print "Did not find $testpath, skipping files in that directory\n";
	} else {
	    foreach my $checkdir (@checkdirs) {
		if (-d "$checkdir") {
		    my $command = 'find ' . $checkdir;
		    my @partial_list = `$command`;
                    foreach my $file (@partial_list) {
                        chomp $file;
			if (-f $file) { push @full_catlist, $file; }
		    }
		}
	    }
	}
    }
    my $keyword_source_file = ''; 
    my $ascii_source = 'catfile_ascii.tmp';
    if (-e "$ascii_source") { unlink($ascii_source); }
    open (my $ascii_fh, '>', "$ascii_source");
    printf "Cataloging %d files in %s\n", scalar @full_catlist, $all_path;
    for (my $i = -1; $i <= $#full_catlist; $i++) {
        my %record;
        # Always put the new CAT file as the first item listed (row 1)
        # Filesize and related characteristics are placeholders that will be revised after
        # the catalog file is completed
        if ($i == -1) {
            %record = (
                'FILENAME', $catdest . '/' . $catfile,
                'FORMAT',   'FITS',
                'TYPE',     'tapecat',
                'FILECLAS', 'catalog',
                'DESCRIP',  "\'FITS-format product catalog\'",
                'FILESIZE', 0,
                'ARCHSIZE', 0,
                'CHECKSUM', 0,
                'GZIP_CRC', 'ffffffff',
                'CKSUM_B4', 0
		);
	} else {
            chomp($full_catlist[$i]);
            if ($full_catlist[$i] =~ /\.cat(\.gz)?$/) { next; } # Skip over cat file in archive
            %record = create_file_characteristics($full_catlist[$i], \@sorted, \%filetype);
            unless (scalar (keys %record)) { next; } # Skip incomplete records/missing files 
            $record{'FILENAME'} = $full_catlist[$i]; 
            $record{'FILENAME'} =~ s/^${all_path}//;
            $record{'FILENAME'} =~ s/\.gz$//;
            # Keyword settings from the first FITS file found
            if ($record{'FORMAT'} eq 'FITS' && $keyword_source_file eq '') { 
                $keyword_source_file = $full_catlist[$i];
	    }
	}
        # write row to FITS file.. [NEED TO MAKE GENERIC BASED ON CDFILE]
        printf $ascii_fh "%64s %16s %32s %32s %64s %d %d %d %8s %d\n",
            $record{'FILENAME'},
            $record{'FORMAT'},
            $record{'TYPE'},
            $record{'FILECLAS'},
            $record{'DESCRIP'},
            $record{'FILESIZE'},
            $record{'ARCHSIZE'},
            $record{'CHECKSUM'},
            $record{'GZIP_CRC'},
            $record{'CKSUM_B4'};
    }
    close $ascii_fh;
    # If no FITS file was found, code cannot get characteristics to populate keywords
    if ($keyword_source_file eq '') {
        #$keyword_source_file = $full_catlist[0];
        print "No FITS for this sequence found, cannot proceed\n";
	exit;
    }
    print "Getting FITS keyword values from $keyword_source_file\n";
    %primary_filename_keywords = keyword_values_fits($keyword_source_file . '[0]', 
                                                 \%primary_filename_keywords);
    %firstext_filename_keywords = keyword_values_fits($keyword_source_file . '[1]', 
                                                  \%firstext_filename_keywords);
    # If $keyword_source_file does not have certain settings, assume time-based settings 
    # from filenames. Some missions do not have all keywords in all files. 
    # NOTE: assumes these are requested keywords: add logic to skip if not requested?
    if (($primary_filename_keywords{'OBS_ID  '} eq '') && ($obsid_date)) { 
        $primary_filename_keywords{'OBS_ID  '} = "'" . $obsid_date . "'";
        $firstext_filename_keywords{'OBS_ID  '} = "'" . $obsid_date . "'";
    }
    if (($primary_filename_keywords{'DATE-OBS'} eq '') && ($date_obs)) { 
        $primary_filename_keywords{'DATE-OBS'} = $date_obs . 'T00:00:00';
        $firstext_filename_keywords{'DATE-OBS'} = $date_obs . 'T00:00:00';
        $primary_filename_keywords{'DATE-END'} = $date_obs . 'T23:59:59';
        $firstext_filename_keywords{'DATE-END'} = $date_obs . 'T23:59:59';
    }
    if (-e "$catfile") { print "Deleting pre-existing local $catfile\n"; unlink($catfile); }
    exec_command('ftcreate cdfile=' . $cdfile . ' datafile=' . $ascii_source . 
                 ' outfile=' . $catfile . ' history=no clobber=yes');
    if ($verbose) { print "ftcreate of $catfile complete, deleting $ascii_source\n"; }
    unlink($ascii_source);
    exec_command('fmodhead ' . $catfile . '+1 ' . $extension);
    exec_command('fmodhead ' . $catfile . '+0 ' . $primary);
    # ftcreate makes TUNIT# fields with ' ' when there are no units: remove these
    exec_command('ftkeypar keyword=TFIELDS infile=' . $catfile . '[1]');
    my $exists = `pget ftkeypar exist`;
    my $column_count = 0; 
    if ($exists =~ /y(es)?/i) {
        my $column_count = `pget ftkeypar value`;
        chomp($column_count);
    }
    for (my $i = 1; $i <= $column_count; $i++) {
        my $tunit = sprintf("TUNIT%d", $i);
        exec_command('ftkeypar keyword=' . $tunit . ' infile=' . $catfile . '[1]');
	$exists = `pget ftkeypar exist`;
	if ($exists =~ /y(es)?/) {
	    my $tunit_setting = `pget ftkeypar value`;
	    chomp($tunit_setting);
	    if ($tunit_setting =~ /^\'\s+\'$/) {
		exec_command("fthedit infile=$catfile operation=delete keyword=$tunit");
	    }
        }
    }
    if ($primary_updated_keywords{'DATE    '}) {
        my $datestring = `date -u +%Y-%m-%dT%H:%M:%S`;
        chomp $datestring;
        $primary_updated_keywords{'DATE    '} = $datestring;
        $firstext_updated_keywords{'DATE    '} = $datestring;
    }
    # Apply the keyword/value pairs to the newcurve file
    assign_keywords($catfile . '[0]', \%primary_fixed_keywords);
    assign_keywords($catfile . '[0]', \%primary_filename_keywords);
    assign_keywords($catfile . '[0]', \%primary_updated_keywords);
    assign_keywords($catfile . '[1]', \%firstext_fixed_keywords);
    assign_keywords($catfile . '[1]', \%firstext_filename_keywords);
    assign_keywords($catfile . '[1]', \%firstext_updated_keywords);
    # Edits to catalog file will have changed its characteristics: revise results
    # Note that the first row in the catalog is always the catalog file itself
    my %record = create_file_characteristics($catfile, \@sorted, \%filetype);
    $record{'FILENAME'} = $catdest . '/' . $catfile;
    foreach my $key (keys %record) {
        exec_command("ftedit infile=${catfile}[1] column=${key} row=1 " .
                     "value=${record{$key}}");
    }
    if ($verbose) { 
        exec_command('ftchecksum ' . $catfile .  '[0] update=yes chatter=2'); 
        exec_command('ftchecksum ' . $catfile .  '[1] update=yes chatter=2'); 
        exec_command('ftverify infile=' . $catfile . ' prstat=yes');
    } else { 
        exec_command('ftchecksum ' . $catfile .  '[0] update=yes chatter=0');  
        exec_command('ftchecksum ' . $catfile .  '[1] update=yes chatter=0');  
        exec_command('ftverify infile=' . $catfile . ' prstat=no errreport=e');
    }
    if (`pget ftverify numerrs` != 0) { die "ftverify found errors in $catfile"; }
    # NOTE: CRC and CKSUM are not updated from compression
    if ($compress_cat) {
        if ($verbose) { exec_command("gzip -v $catfile"); }
        else          { exec_command("gzip $catfile"); }
        $catfile .= '.gz';
    }
    if ($archive_cat) {
        if ($verbose) { print "Moving $catfile to archive\n"; }
        exec_command("mv ${catfile} ${all_path}${catdest}/${catfile}");
    }
}

 
# Replace all variable names with patterns for pattern file/directory searchs
sub replace_varname_with_pattern {
    my $search_string = shift;
    my ($pattern_definitions_ref) = @_;
    my %pattern_definitions = %{$pattern_definitions_ref};

    my @fields = split(/:/, $search_string);
    my $num_variables = int(scalar (@fields) - 1) / 2;
    my @vars;
    my @patterns;
    my $pattern_string = $search_string;
    for (my $i = 0; $i < $num_variables; $i++) {
        my $index = $i * 2;
        my $variable_name = $fields[$index+1];
        $vars[$i] = $variable_name;
        if ($pattern_definitions{$variable_name}) { 
            my $replacement_pattern = $pattern_definitions{$variable_name};
            $pattern_string =~ s/\:${variable_name}\:/${replacement_pattern}/; 
            $patterns[$i] = $replacement_pattern;
        } else {
            die "No search pattern found for variable $variable_name"; 
        }
    }
    return($pattern_string, \@vars, \@patterns);
}

# Replace all variable names with patterns for pattern file/directory searchs
sub find_value_with_pattern {
    my $search_string = shift;
    my $resolved_string = shift;
    my ($vars_ref,$patterns_ref) = @_;
    my @vars = @{$vars_ref};
    my @patterns = @{$patterns_ref};
    my %var_setting;

    my @fields = split(/:/, $search_string);
    my $num_variables = int(scalar (@fields) - 1) / 2;
    my $pattern_string = $search_string;
    my $partial = $fields[0];
    for (my $i = 0; $i < $num_variables; $i++) {
        my $index = ($i * 2) + 1;
        my $variable_name = $fields[$index];
	my $pattern;
        for (my $j = 0; $j <= $#vars; $j++) {
            if ($vars[$j] eq $variable_name) {
                $pattern = $patterns[$j];
                last;
            }
        }
        if ($resolved_string =~ /^${partial}(${pattern})/) {
            $var_setting{$variable_name} = $1;
            $partial .= $1;
            if ($fields[$index+1]) { $partial .= $fields[$index+1]; } 
        } else { die "Could not find variable"; }
    } 
    return(%var_setting);
}


#########################################################
# Subroutine to concatenate date and time strings:
# Input:
# - date_in (Str): Date component of datetime (YYYY-mm-dd)
# - time_in (Str): Time component of datetime (HH:MM:SS) 
#
# Output:
# - dt_out (Str):  Full datetime output (YYYY-mm-dd HH:MM:SS)
#########################################################
sub dt_concatenate {
    my ($date_in, $time_in) = @_;
    my $dt_out = "${date_in} ${time_in}";
    return $dt_out;
}


#########################################################
# Subroutine to concatenate date and time strings in ISO format
# Input:
# - date_in (Str): Date component of datetime (YYYY-mm-dd)
# - time_in (Str): Time component of datetime (HH:MM:SS) 
#
# Output:
# - dt_out (Str):  Full datetime output (YYYY-mm-ddTHH:MM:SS)
#########################################################
sub iso_concatenate {
    my ($date_in, $time_in) = @_;
    my $dt_out = "${date_in}T${time_in}";
    return $dt_out;
}


# Get the values for the requested keywords from the file $fits_file

sub keyword_values_fits {
    my ($fits_file, $keyword_ref) = @_;
    my %keywords = %{$keyword_ref};

    foreach my $key (keys %keywords) {
        exec_command('ftkeypar ' . $fits_file . ' ' . $key);
        my $value = `pget ftkeypar value`;
        chomp($value);
        $keywords{$key} = $value;
	my $exists = `pget ftkeypar exist`;
	if ($exists =~ /no/i) { $keywords{$key} = undef; }
    }
    return (%keywords);
}

# Get the catalog file characteristics (checksums, file size, etc) for file $infile
sub create_file_characteristics {
    my ($infile, $sorted_ref, $filetype_ref) = @_;
    my @sorted_types = @{$sorted_ref};
    my %filetype = %{$filetype_ref};
    my %results;
	
    my $cksum    = undef;
    my $size     = undef;;
    my $cksum_b4 = undef;
    my $filesize = undef;;
    my $line = `cksum $infile `;
    my @results = split /\s+/, $line;
    $cksum = $results[0];
    $size  = $results[1] / 1024;
    if ( !defined $cksum || !defined $size ) {
        print "No output from cksum of $infile\n";
    }
    # Get filesize before gzip by piping gunzip (gzip -d) to STDOUT (-c)
    # WRITING TO DISK WOULD GUNZIP FILES IN THE ARCHIVE (which would be bad)
    if ($infile =~ /.gz$/) {
        $line = `gzip -d -c $infile | cksum  `;
        @results = split /\s+/, $line;
        $cksum_b4 = $results[0];
        $filesize  = $results[1] / 1024;
    } else {
        $cksum_b4 = $cksum;
        $filesize = $size;
    }
    if ( !defined $cksum_b4 || !defined $filesize ) {
        print "No output from \"gunzip\" cksum of $infile\n";
    }
    my $crc;
    my $gzip_size;
    if ($infile =~ /.gz$/) {
        open GZIP, "gzip --list --verbose $infile |";    
        if ( $?) { print "$infile not understood by gzip\n"; }
        <GZIP>; 
        $_ = <GZIP>;
        my @line = split /\s+/, $_;
        $crc = $line[1];
        $gzip_size = $line[5] / 1024;
        $filesize = $line[6] / 1024;
        close GZIP;
    } else {
        $crc = 'ffffffff'; # CRC is "-1" if file is not GZIPed
        $gzip_size = $size;
        $filesize = $size;
    }
    my $filetest = $infile;
    $filetest =~ s/\.gz$//;
    my $found = 0;
    for (my $i = 0; $i <= $#sorted_types; $i++) {
        my $key = $sorted_types[$i];
        if ($filetest =~ /${key}$/) {
            $results{'FORMAT'}   = $filetype{$key}[0];
            $results{'TYPE'}     = $filetype{$key}[1];
            $results{'FILECLAS'} = $filetype{$key}[2];
	    $results{'DESCRIP'}  = $filetype{$key}[3];
            $found = 1;
            last;
	}
    }
    $results{'FILESIZE'} = ceil($filesize); # Round file size up with POSIX ceiling
    $results{'ARCHSIZE'} = ceil($gzip_size);
    $results{'CHECKSUM'} = $cksum;
    $results{'GZIP_CRC'} = $crc;
    $results{'CKSUM_B4'} = $cksum_b4;
    unless ($found) { 
        print "Failed to find $filetest file type, skipping\n";
        %results = ();
    }
    return %results;
}

# Subroutine to apply keyword values already set in an associative array
# to a newly created FITS file which has placeholder values for all these keyword
# pairs.

sub assign_keywords {
    my ($filename, $keyword_ref) = @_;
    my %keywords = %{$keyword_ref};

    # Add ' to force literal interpretation of filename if FITS filename includes a 
    # [] extension indicator
    if ($filename =~ /\[\d+\]$/) { $filename = "\'${filename}\'"; }
    foreach my $key (sort keys %keywords) {
        my $value = $keywords{$key};
	if (defined $value) {
            # Some values look like numbers, but need to be strings (e.g. OBS_ID)
            if ($value =~ /\'\d+\'/) {
                exec_command("fparkey value=\"" . $value . "\" fitsfile=" . 
                         $filename . ' keyword=' . $key); 
	    } else {
                exec_command("fparkey value=" . $value . " fitsfile=" . 
                         $filename . ' keyword=' . $key); 
	    }
	} else {
	    print "WARNING! Expected a value for keyword $key in $filename, but not defined!\n";
	}
    }

}

# Execute system commands, collecting error status and printing to STDOUT if
# verbose is enabled, and printing out an error alert is something went awry

sub exec_command {
    my ($instruction) = @_;
    
    if ($verbose) { print "$instruction\n"; }
    my $status = system($instruction);
    if ($status != 0) { 
        print "Error with $instruction\n"; 
        if ($verbose) { die "Exiting..."; }
    }
    return($status);
}    

# Subroutine to read a template described the keyword to be placed in the final
# FITS file. The template should have three sections: keywords that are fixed, 
# keywords that will updated by contents of the source files, and keywords that
# be updated at the end of processing (usually processing date (now: DATE),
# DATASUM and CHECKSUM). The sections are seperated from each other by comment
# lines (Lines that start with "#", blank lines, or (most often) both. 
# Lines in the template file should look much like the output of fdump. E.g.
# "TELESCOP= 'MAXI    '           / Mission name"
# with FITS keyword, '=', value, '/', and a comment string. The comment string
# is ignored at this stage.

sub read_keyword_template {
    my $template_file = shift;
    my %fixed_set;
    my %toalter_set;
    my %toupdate_set;
    my $dataline = 0;
    my $stage = 0;

    open(my $fh, '<', $template_file);
    while (my $line = <$fh>) {
        if (($line =~ /^\s+$/) || ($line =~ /^#/)) {
            $dataline = 0;
            next;
        }
        if ($line =~ /\=/) {
            # If preceeding line(s) are/were either blank or comments, shift to next 
            # set of keywords on first new line with keyword data
            if ($dataline == 0) { $dataline = 1; $stage++; } 
            my ($keyword_value,$comment) = split(/ \/ /,$line);
            my ($keyword, $value) = split(/\=/, $keyword_value);
            $value =~ s/^\s+//;
            $value =~ s/\s+$//;
            if ($stage == 1)    { $fixed_set{$keyword} = $value; }
            elsif ($stage == 2) { $toalter_set{$keyword} = $value; }
            elsif ($stage == 3) { $toupdate_set{$keyword} = $value; }
            else { die "Found more than three sections of keywords in $template_file"; }
	}
    }
    close $fh;
    return(\%fixed_set, \%toalter_set, \%toupdate_set);
}

# Standard Usage message

sub use_message {
    print <<EndOfMessage
Version: $VERSION

Usage: make_generic_cat_file.pl --template cat_template --primary primary_template \
                                --extension 1stext_template [--cdfile cdname] \
                                [--range VARIABLE --min MINIMUM --max MAXIMUM]
                                [--verbose] [--archive] [--help]
WHERE
   - cat_template       is a template file for what types of files are in the archive
                        and should be listed in the CAT file, and their characteristics
   - primary_template   is a template file listing the FITS keywords to go into the 
                        primary array of the catalog file. This should have three sections
                        seperated by either a blank line or comment line (or both)
                        for keywords that will be fixed in all catalogs (e.g. TELESCOP), 
                        keywords whose values are read from source FITS files in the 
                        archive (e.g. DATE-OBS), and keywords that will be updated by
                        this script (CHECKSUM, DATASUM, and DATE (processing date)
  - 1stext_template     is a template file listing the FITS keywords to go into the 
                        first binary extension where the catalog will be listed. It should
                        be organized in the same way as the primary array template
  - cdname              is an optional list of the columns and formats for the catalog
                        listing. 
  --range               allows a single variable (e.g. OBSID) to be limited to a specified
                        range of values. Overrides any range settings in the cat_template.
                        MUST include --min and --max if --range is invoked.
  - archive             is an optional flag to specify that the CAT files made by the 
                        script should be moved directly into the archive once created:
                        default is to write to the local directory.
  - verbose             is an optional flag to run this script in a more chatty mode
                        with more screen output, useful in debugging problems.
  - help                Prints this help message 
EndOfMessage
}
