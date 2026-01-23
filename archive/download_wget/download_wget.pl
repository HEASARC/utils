#!/usr/bin/perl
#
# Usage:    download_wget.pl URL [debug]
#
#           Where URL is a URL to the HEASARC FTP area 
#           in the URL e.g., https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/ 
#
#  Version 1.0    J. Allen 2019  (LA/TM)
#          1.001  Changed 'index' to 'index*' to eliminate all index.html creation (2020/06/26)
#                 Check if any files are downloaded: report if nothing found (possible URL error)

use strict;
use warnings;
my $command;
my $status;

# Give the help is the script is invoke with no argument or too many 
if (($#ARGV < 0) || ($#ARGV > 1)) { usage_message(); }

# Given help if "download help" 
if ($ARGV[0] =~ /help/i) { usage_message(); }

# Exit with usage message if given a non-HEASARC URL
if ($ARGV[0] !~ /^https:\/\/heasarc\.gsfc\.nasa\.gov/) { 
    print "ERROR: This is not a HEASARC server address\n\n";
    usage_message(); 
} 

# Assign first command argument to $url
my $url = $ARGV[0];

# Hidden debug feature invoked as "download_wget.pl url debug". Exit with help for anything else
# Debug is shared with subroutines, so "our" instead of local "my"
our $debug = 0;
if ($#ARGV > 0) {
    if ($ARGV[1] =~ /debug/i) { $debug = 1; }
    else { usage_message(); }
}                                                    

# In debug mode, print the full URL provided
if ($debug) { print " Version 1.0001 (26 Jun 2020 version)  \n" ; } 
if ($debug) { print " 1- $url \n" ; } 

# Find the part of the URL after the protocol.
my $xurl = $url; 
$xurl =~ s,(^https?|ftp)://,,;

# Look how many  /
my @flds = split("/", $xurl);
if ($debug) { print " 3- $#flds \n"; } 
if ($debug) { print " 4- $flds[$#flds] \n"; }
my $cut_count = $#flds - 1; 
if ($debug) { print " 5- cut-dir=$cut_count \n"; }
if ($#flds < 3) { print "URL $url will try to download too much of the archive, exiting."; exit; }
elsif ($#flds == 3) { 
    print "WARNING! Valid URL, but you may be getting too much data...\n"; 
    print "Please try again with a more specific pattern match request\n";
    exit;
}

# Look for certain wildcards in URL ([(stuff)], ?, and *)
# Call subroutine to deal with wildcards else proceed
if ($url =~ /\*|\[.+\]|\?/) { wget_wildcard_search($url); }

# Add a '/' to the end of the URL to test if this is directory, modify URL if directory
unless ($url =~ /\/$/) {
    my $testurl = $url . '/';
    # Silent URL test (--spider tests validity of URL, but does not download)
    $command = "wget --spider --quiet $testurl"; 
    $status = system($command);
    if ($status == 0) {
        # It is a directory: change the URL to include the '/'
        $url = $testurl;
        if ($debug) { print " 2- $url \n";} 
    }
}
# Run the download command (-q/--quiet unless debug is set)
$command = "wget -q -nH -r -c -N --retr-symlinks -e robots=off -N -np -R \'index*\' " .
           "--cut-dirs=$cut_count $url";
if ($debug) { 
   $command =~ s/^wget \-q/wget/;
    print "$command\n"; 
} else {
    print "Downloading $url\n"; 
} 
$status = system($command);
wget_check_status($status);
if ($status == 0) {
  print "Download complete.\n";
  my $local_dir = $flds[$#flds];
  # Check up to three directories down in download(s) for files
  my @results = glob "$local_dir/*";
  my @sub_results = glob "$local_dir/*/*";
  push @results, @sub_results;
  if ($#sub_results >= 0) {
    @sub_results = glob "$local_dir/*/*/*";
    push @results,@sub_results;
  }
  my @files;
  foreach my $found (@results) {
    if (-f "$found") { push @files, $found; }
  }
  if ($#files < 0) {
    print "Download found no files in $local_dir: were there no files in $url?\n";
  }
}
else { print "Download failed: please check URL\n"; }

#====== END MAIN ROUTINE ========

sub wget_check_status {
    my ($status) = @_;

    if ($status == -1) {
        print "$command failed to execute\n";
    }
    elsif ($status & 127) {
        printf "wget died with signal %d, %s coredump\n",
            ($status & 127),  ($status & 128) ? 'with' : 'without';
    }
}

sub get_files_from_index {
    my ($indexfile) = @_;
    my @results = ();

    if (-e "$indexfile") {
        open (my $fh, '<', "$indexfile");
        while (<$fh>) {
            unless (/^<img src/) { next; }
            chomp;
            my $file;
            if (/a href/) {
                my @elements = split('a href');
		my $start_index = index($elements[1],"=\"");
		my $stop_index  = index($elements[1],"\">");
                $file = substr($elements[1],$start_index+2,$stop_index-($start_index+2));
	    }
            push (@results, $file);
	}
        close $fh;
    }
    return \@results;
}

sub wget_wildcard_search {
    my ($url) = @_;

    my $status;
    my $command;
    my $recurse_limit = 5; # Look up to 5 directories down from URL in search
    my @files = ();        # List where all the files matching the pattern will go

    my $testurl = $url;
    $testurl =~ s,^(https?|ftp)://,,;
    my $headstring = $1;
    # Split URL into directories and find the first one with a wildcard
    # $first_wild is the index of this first wildcard
    # $wild_count is the number of wildcards in the URLs: for multiple wildcards,
    # all files will be downloaded to the current directory
    # E.g. https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2019_01/*/auxil/*.mkf.gz
    # would put all the .mkf.gz files into the current directory
    # If there is only a single wildcard setting, maintain the directory structure
    # E.g. https://heasarc.gsfc.nasa.gov/FTP/maxi/data/obs/MJD58000/MJD5854[2-4]
    # would create MJD58542, MJD58543, and MJD58544 and populate those with all the 
    # matching subdirectories
    my @fields = split('/',$testurl);
    my $first_wild = $#fields-1;
    my $wild_count = 0;
    for (my $i = $#fields; $i >= 0; $i--) {
        if ($fields[$i] =~ /\*|\[.+\]|\?/) { $first_wild = $i; $wild_count++; }
    }
    if ($first_wild <= 3) { print "Wildcard in URL $url too close to root, not allowed\n"; exit; }
   
    my $current_wild = 999;
    my $previous_wild = 999;
    for (my $i = 0; $i <= $#fields; $i++) {
        if ($fields[$i] =~ /\*|\[.+\]|\?/) { $current_wild = $i; }
        if ($previous_wild == ($current_wild - 1)) {
            print "URL $url has two consecutive wildcards: not allowed (too much data)"; 
            exit;
	} else { $previous_wild = $current_wild; }

    }

    # Convert unix-style wildcards into Perl regex for pattern matching
    # $search is the full pattern, @search is the pattern within each level
    # of the directories.
    # E.g. https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_09/*/auxil/ni*.att.gz
    # $search = '.*/auxil/ni*.att.gz'
    # @search_fields = ('.*','auxil','ni.*.att.gz')
    # @search_fields allows pattern matching within each directory
    my $search = '';
    my @search_fields;
    for (my $i = 0; $i <= ($#fields-$first_wild); $i++) {
        my $temp = $fields[$i + $first_wild];
        $temp =~ s/\?/\./;   # Convert ? to . for Perl regex
        $temp =~ s/\*/\.\*/; # Convert * to .* for Perl regex
        push(@search_fields,$temp);
        $search .= $temp . '/';
    }
    chop($search); 

    # Complete regex style search string
    if ($debug) { print " search string - $search \n" ;} 

    # Build a list of all files and directories that match the pattern
    # First, make the base URL
    my $tmpurl = $headstring . '://'; 
    for (my $i = 0; $i < $first_wild; $i++) { $tmpurl .= $fields[$i] . '/'; }
    chop($tmpurl);

    my @search = ('/');
    my @next_search = ();
    print "Building a list of all files and directories that match the pattern search\n";
    my $searching = 1;
    while ($searching) {  
        my $dir_count = 0;
        foreach my $dir (@search) {
            # List all files and directories in $dir in index.html, and 
            # sort these into directories that match pattern and files that match 
            # pattern and ignore all others. Directories that match pattern are 
            # stored in @next_search and will be traversed to look for files in
            # the next step down the directory structure ($i+1). Continue until
            # running out of files or hit the recursion limit (limit on how far down the 
            # directory structure to look)
            my $testurl = $tmpurl . $dir;
            if ($debug) { print "Searching $testurl\n"; }
            elsif ((($dir_count % 20) == 0) && ($dir_count != 0) && (scalar(@search) > 20))
	        { printf "%d of %d directories searched\n", $dir_count, scalar(@search); }
            my @search_depth = split('/',$dir);
            my $current_level = $#search_depth;
            if ($current_level < 0) { $current_level = 0; }
            $command = 'wget --no-parent --no-host-directories ' .
                       "--quiet --retr-symlinks -e robots=off " .
                       "--no-remove-listing $testurl";
            if (-e "index.html") { unlink "index.html"; }
            $status = system($command);
            my $list = get_files_from_index("index.html");
            unlink "index.html";
            $dir_count++;
            # Check which items are files or directories, check if they match the 
            # search pattern, and add any matching files to @files, and 
            # any matching directories to @next_search to set up a search 
            # (e.g. url is 
            # https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_09/*/auxil/ni*.att.gz
            # all directories below 2018_09 will match, need to be searched for auxil,
            # and all files in auxil matching ni*.att.gz listed to @files)
            my $regex = '.*';
	    if ($current_level <= $#search_fields) { $regex = $search_fields[$current_level]; }
            my $next_dir = '';
            my @new_files = ();
            my $next_level = $current_level + 1;
            if ($next_level <= $#search_fields) { 
                $next_dir = $search_fields[$next_level] . '/'; 
                # If the next directory has no wildcards, we can skip ahead in search
                # but if it has any special characters, perform a search through the 
                # directory tree
                if ($next_dir !~ /^\w+\/$/) { $next_dir = ''; } 
            }
            if (scalar (@$list == 0)) { print "Download failed: check URL\n"; exit; }
            foreach my $listing (@$list) {
                if ($listing =~ /\/$/) {
                    # This is a directory: If match, push to next search
                    if ($listing =~ /${regex}/) {
                        push(@next_search, $dir . $listing . $next_dir);
		    }
	        } else {
                    # This is a file: record the full URL if it matches the pattern
                    if ($listing =~ /${regex}/) { push(@new_files, $testurl . $listing); }
   	        }
	    }
            if (scalar(@new_files)) {
                if ($debug) {
                    if (scalar(@new_files) == 1) 
                        { printf "Found %d file\n", scalar(@new_files); }
                    else
                        { printf "Found %d files\n", scalar(@new_files); }
		}
	    }
            push(@files,@new_files);
	} # End loop through current directory: hand @next_search to @search, or end
          # search if no more directories to search were found
        if (scalar(@files)) { 
            if (scalar(@files) == 1) 
                { printf "Found %d file to download so far.\n", scalar(@files); }
            else 
                { printf "Found %d files to download so far.\n", scalar(@files); }
	}
        @search = ();
        if ($#next_search != -1) { 
            @search = @next_search; 
            @next_search = (); 
            printf "Found %d directories\n", scalar(@search);
            if (scalar(@search) > 100) { print "Please be patient: this may take a while\n"; }
        } else { $searching = 0; } # Terminate loop when no further matches are found
    } # End while searching loop
        
    print "Found all matching files, now downloading...\n";
    my $cnt = 0;
    if ($#files < 0) { print "No files to download: check url\n"; exit; }
    foreach my $download (@files) {
        # Run WGET on all files that match pattern, then exit
        $command = "wget -q -nH -r -c -N --retr-symlinks -e robots=off -N -np -nd -R \'index\' " .
               sprintf("--cut-dirs=%d", $first_wild-1) . " $download";
        if ($debug) { $command =~ s/^wget \-q/wget/; print "$command\n"; }
        else {
            if (($cnt % 20) == 0) { printf "%d of %d files downloaded.\n", $cnt, scalar(@files); }
	}
        $cnt++;
        if ($wild_count == 1) { $command =~ s/ -nd//; } # Retain directories if only one wildcard
        $status = system($command);
        wget_check_status($status);
        unless ($status == 0) { print "Download failed: check that $download is a valid address\n"; }
    }
    if (($cnt % 20) != 0) { printf "%d of %d downloaded.\n", $cnt, scalar(@files); }
    print "Downloads completed.\n";
    exit;
} # end wild card search 
 

sub usage_message {

    print << "EndOfHelp";

Usage :  download_wget url 
 
Downloads data to the local computer from the HEASARC archive using wget command and a url. 
The possible type of downloads are :
 a) single directory corresponding to a specific observation
 b) range of directories corresponding to a range of observations 
 c) single file 
 d) type of files 
 
a) The command for a single directory is : 

 >  download_wget.pl https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/

  On the local computer the downloaded data are in the directory 1050020180/ maintaining  
  with the archive structure.

b) The command for a range of directories is :

 > download_wget.pl "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/105002018[01]"

  On the local computer the downloaded data are in the directories 1050020180/ and 
  1050020181/. The archive structure is maintained within each directory.

c) The command for a single file is :

 > download_wget.pl https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/auxil/ni1050020180.att.gz

  On the local computer the file is downloaded in the directory where the script is invoked. 

d) The command for multiple files with the same pattern : 

 > download_wget.pl "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/*/auxil/ni*.att.gz"

  On the local computer the files are downloaded in the directory where the script is invoked.
  The archive structure is not maintained.


Within a given url wildcard are allowed to search specific files that are located in different 
directories with a well defined pattern. The wildcard completion allowed are * and [ ]. A 
maximum of two non consegutive wildcards are allowed for url.

If the transfer is interruped, the user may restart the same command in the same directory and only
the remaining files are downloaded.  

EndOfHelp
    exit;
}

#
# 
#Allowed                                                                                                                                  
# >  download_wget https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/
# >  download_wget "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/105002018[01]"
# >  download_wget "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/*/auxil/ni*.att.gz"
#                               
# Not allowed and need a more specific pattern
# > download_wget "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/*/*/ni*.att.gz"
# > download_wget "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/
