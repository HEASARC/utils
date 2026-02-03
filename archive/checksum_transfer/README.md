# checksum transfer

Instructions for using the checksum transfer script to
compute checksums for a set of files at their origin, and to
test that the checksum matches on arrival at the destination

## Contents

`checksum_transfer.pl`
Perl script to generate (at origin and destination) checksums,
and (at destination) compare and checksums

`checksum_origin.txt`
Sample list of files and checksum values generated at origin
(example from Halosat; origin at U. Iowa 2020-02-07)

`checksum_destination.txt`   
Sample list of matching destination files and checksum values
generated at destination (example from Halosat; destination
at HEASARC 2020-02-07)

`dirlist.txt`
Example ascii file listing directories

## PURPOSE

Verify that a data set has not been truncated accidentally during transfer, by 
comparing starting and finishing checksums.

## CALLING SEQUENCE

```
checksum_transfer.pl [TOPDIR FILE] [LOCATION] [ORIGIN FILE] [DESTINATION FILE]
```

## INPUT

These parameters will be prompted for if not provided.

TOPDIR FILE: ascii file listing full paths to all top directories to be checked.

LOCATION: location at which checksum is being done; "O" for origin, "D" for destination

ORIGIN FILE: name of file containing checksums from origin location (Origin mode: will be 
written, Destination mode: will be read)

DESTINATION FILE: name of file containing checksums from destination location (only used 
in Destination mode)


## OUTPUT

an ascii text file listing recursively every file under the input top directory list and 
its checksum values.  If doing a destination check, this file is compared against the 


## Example

At the origin

```
% cd /path/to/local/files
% find . -maxdepth 1 -type d > dirlist.txt
% checksum_tranfer.pl dirlist.txt O checksum_origin.txt
```

output is checksum_origin.txt

The operator at origin delivers the files and the dirlist.txt and checksum_origin.txt file to the
destination operator.


At the destination, the operator places the files in their destination directory, then

```
% cd /path/to/destination/files
% checksum_tranfer.pl dirlist.txt D checksum_origin.txt checksum_destination.txt
```

The output is checksum_destination.txt. The code will report

Files checksum_origin.txt and checksum_destination are identical; all checksums verified.

or will give an listing of all differences found.


## Author

Jesse Allen (NASA/GSFC HEASARC)
18 September 2024

