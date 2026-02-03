# Creating a Generic Catalog File

make_generic_catfile.pl is a Perl script to make "any mission" archive 
catalog files.          

## HEASoft

This tool requires a working FTOOLS set up. The command 

```
printenv HEADAS
```

returns the HEASoft installation directory: if nothing is returned,
install or enable FTOOLS before proceeding further.


## Command

To run the catalog generation script, type:

```
./make_generic_catfile.pl --cdfile cdfile_cat_basic.txt \
                          --primary template_primary.txt \
                          --extension template_1stext.txt \
                          --template catfile_template.txt
                          --archive --compress
```

where

- `cdfile_cat_basic.txt` is a standard column description ASCII list used to generate the file listings in the CATALOG binary table extension in the .cat file
- `template_primary.txt` is a standard template for FITS keywords, values, and comments to be placed in the primary array header of  the .cat file      
- `template_1stext.txt` is a standard template for FITS keywords, values, and comments to be placed in the CATALOG binary table extension of the .cat file
- `catfile_template.txt` is a file that describes all the different types of files to be found for each observation, plus search pattern settings for finding them

The code makes a series of .cat FITS format catalog files by looking in 
the data archives for all files that match patterns that are specified in the 
`catfile_template.txt` file. Each CAT file should have, at a minimum, columns 
(described in the column description `cdfile.txt` file), include all file 
names, file sizes (uncompressed and, if they are provided in gzip form in the 
archive, compressed), gzip CRC checksums, checksums (based on uncompressed 
and compressed sizes), a basic short file description, a file format (usually 
FITS), type, and class. In instances where the file is stored in the archive 
uncompressed, the gzip CRC is set to `ffffffff` (-1 in hexidecimal) and the
compressed checksum and filesize are set to be the same as the uncompressed. 
The CAT file should also include an entry for itself, based on filesizes and 
checksums once it is fully created.

In addition to the required templates, there are settings to control the 
behavior of the script

`--archive`              will write the .cat file into the observation data 
                       directories in the subdirectory set by the variable 
                       CATDEST (or CATDIR) in the catfile_template.txt. If not
                       set, the .cat file will be written in the current directory

`--compress`             will compress the .cat file (using gzip) 

`--help`                 will give a summary of all possible required and optional fields
                       to the script

`--verbose`              will generate a more "talkative" output that can be useful 
                       should problems arises

`--range variable MIN MAX` allows one variable to be restricted to a smaller range
                      of possible values. This can be useful when generating cat
                      files for part of an archive (e.g. an active archive where 
                      cat files are needed only for new date)
                      E.g. if there is a variable OBSID in the catfile_template.txt
                      of the form YYMMDD, "--range OBS_ID 220101 220110"
                      would restrict the script to only generating new .cat files
                      between these values.


## FORMAT OF THE CATALOG FILE: cdfile_cat_basic.txt

This is an ascii file containing a list of column names that will be included in the
cat file and their format. Example :

```
FILENAME 64A
FORMAT   16A
TYPE     32A
FILECLAS 32A
DESCRIP  64A
FILESIZE 1V kilobytes
ARCHSIZE 1V kilobytes
CHECKSUM 1V
GZIP_CRC 8A
CKSUM_B4 1V
```

This is identical to the format used for column description files used by
the FTOOL ftcreate.

These columns are:

`FILENAME` 
a column which contains the names of files, including all relative paths
within the individual observation (e.g. auxil/file.fits)

`FORMAT`
A single word description of the format of the file. FITS is most common, but others 
such as HTML or GIF may arise in some archives

`TYPE`
A short word describing the file type, usually a truncated description (e.g. "event", 
"orbit", "ancillary", or "housekeeping")

`FILECLAS` 
A file classification word. Often this is the same as TYPE, though in some instances 
TYPE is more specific (e.g. FILECLAS might be "housekeeping" while TYPE is "det1_hk" for 
an individual detector housekeeping file).

`DESCRIP`
A short one-line file description 


## FORMAT FOR THE PRIMARY: template_primary.txt

The primary array template includes three sections of FITS keywords to be 
placed in the primary array header for the catalog file. The first set is 
fixed keywords (values are the same in all catalog files), the second set
(marked off from the first by a comment line) shows values that will 
change from one catalog to the next and are based on values for the FITS 
keywords in files for each observations, and the third set are data quality
tracking settings set by the catfile creation. These last two sets are shown
here as blank settings

```
# Template for BurstCube FITS catalog file primary header
# 2021-10-28
#
# standard valid for all files
TELESCOP= 'BURSTCUBE'          / Mission name
INSTRUME= 'CS       '          / Instrument name 
ORIGIN  = 'GSFC    '           / Origin of FITS file

# Placeholder settings, computed from data file(s)
OBS_ID  = ''                   / Observation ID
DATE-OBS= ''                   / Start date of observations
DATE-END= ''                   / End date of observations

##Recalculated when .cat file created
DATE    = ''                   / date of file creation (GMT)
CHECKSUM= ''                   / HDU checksum 
```

## FORMATS FOR THE 1st EXTENSION: template_1stext.txt

The first binary table extension of the catalog file contains the catalog information.
Typically, this is the same as the primary array, but with EXTNAME and DATASUM
added:

```
# Template for BurstCube FITS catalog file primary header
# 2021-10-28
#
# standard valid for all files
TELESCOP= 'BURSTCUBE'          / Mission name
INSTRUME= 'CS '                / Instrument name
ORIGIN  = 'GSFC    '           / Origin of FITS file
EXTNAME = 'CATALOG '           / name of this binary table extension

# Placeholder settings, computed from data file(s)
OBS_ID  = ''                   / Observation ID
DATE-OBS= ''                   / Start date of observations
DATE-END= ''                   / End date of observations

##Recalculated when .cat file created
DATE    = ''                   / date of file creation (GMT)
CHECKSUM= ''                   / HDU checksum 
DATASUM = ''                   / data unit checksum 
```

Additional FITS keywords can be added to either or both of the primary and 
first extension of the catalog file if desired.


## FORMAT FOR THE CATALOG file: catfile_template.txt

It is an ascii files containing the list of files and their attributes
as well as a set of keywords that defines file-pattern search to list
all files in an observation and categorize them.
The file contains three parts. The first part set pattern via the 'variables' token
to describe the directory path and if necessary any file pattern that can not 
be recontructed with the path. The second part set search for file via pattern 
search set in 'filelist' token. The last part set the single file attributes.

### First part: setting of the "variables" token

In a typical archive, these is a top level path which is the same for all
observations (PATH or ROOTPATH), and a subdirectory below which varies between
observations (SUBDIR). Often this is a multi-level structure (e.g. 
"2020_01/689900/" for an observation "689900" obtained in January 2020)
Each observation may be in a single directory or a nested structure with 
multiple directories (e.g. "events", "images", "auxil"). 

The catfile_template requires a setting for PATH, SUBDIR, and CATNAME.
However, it is often useful for a mission-specified variable or multiple
variables to be defined. These use Unix pattern matching (not regex or 
egrep patterns) to match directories and filenames.

### Second part: setting of the "filelist" token 

The list of files is defined by the setting "filelist" which uses Unix
pattern matching (NOT regex or egrep!) to generate a list of all files
for an observation.

The resulting list of files is matched by the "category" setting using 
file suffixes. Multiple files can match a single category, but the suffix
must be able to identify all file types present: if, for example, there
are files "688900.hk", "688900_det1.hk", and "688900_det2.hk", 
using "*.hk" would match all three. Specifying "_det1.hk", "_det2.hk" and
".hk" would distinguish between them.

### Third part: setting the file attributes

The third part contains the list of the file attributes. These are file suffix 
(the last character after the dot in the filename), FORMAT (type of format as 
for FITS, ASCII. etc.), type (describe what kind of file), file class 
(define the file class) , and description (what the file contains).

The catfile_template can contain comment lines, indicated by a line
starting with "#". Blank lines are ignored.

Example: For mission ZZEXPLORER, data are placed into the archive in 
subdirectories of the format YYYY_MM/XXXXXX where YYYY is the year, MM is the month, 
and XXXXXX is a 6-digit unique sequence number for each observation. Each observation 
contains a number of files, sorted into different subdirectories within 
each observation. The catfile, once generated, will be placed in the 
"auxil" directory. There is also multiple possible subdirectories in the event/ folder
with their own ID pattern for the directories. ZZEXPLORER data has dates no earlier
than January 2000 and no later than December 2099.

### Sample catfile template file

```
variable PATH /path/to/archive/
variable OBS_ID [0-9][0-9][0-9][0-9][0-9]
variable SECOND_ID [0-9][0-9][0-9][0-9][0-9][a-z]
variable SUBDIR 20[0-9][0-9]_[01][0-9]/:OBS_ID:/

# Additional mission specific settings, if any 
# ZZEXPLORER has four detectors with names a, b, c, and d
variable DETNAME [a-d] 

# Describe the catfile name and its destination 
variable CATNAME zz:OBS_ID:.cat
variable CATDEST auxil

# Describe the filename pattern for all files in a single observation
filelist auxil/zz:OBS_ID:*
filelist events/zz:OBS_ID:*
filelist events/:SECOND_ID:/zz:OBS_ID:*
filelist spectra/zz:OBS_ID::DETNAME:.pha
filelist log/zz:OBS_ID:*.html
#
# Template for file types, file class, and description
#
#        file_suffix       FORMAT TYPE         FILECLASS     DESCRIPTION
#
category .cat               FITS  tapecat      catalog       FITS-format product catalog
category .att               FITS  attitude     attitude      Attitude file
category .orb               FITS  orbit        orbit         Orbit file
category .mkf               FITS  filter       filter        Filter file
category _uf.evt            FITS  uf_evt       event         unfiltered event file
category a.pha              FITS  zza_pha      spectra       Detector A PHA spectra file
category b.pha              FITS  zzb_pha      spectra       Detector B PHA spectra file
category c.pha              FITS  zzc_pha      spectra       Detector C PHA spectra file
category d.pha              FITS  zzd_pha      spectra       Detector D PHA spectra file
category _joblog.html       HTML  joblog       log           HTML processing log file
category _errlog.html       HTML  errlog       log           HTML error log file
```

The final catfile generated by this code has one row per file for the given 
observation, with all the file information as defined by the column decription 
file (FILECLAS, FILENAME, checksums, etc).

If any files are found that do not match any of the possible suffix patterns, 
make_generic_catfile.pl prints a warning message alerting the operator: this 
may indicate unintended temporary files not cleaned out before archiving, or
incomplete file suffix lists.


## How to use the generic catalog file generator for CALET CGBM

This describes the steps for creating .cat files for the CALET CGBM data

0/ It requires a working FTOOLS set up. The command 

```
printenv HEADAS
```

returns the HEASOFT installation directory: if nothing is returned,
install or enable FTOOLS before proceeding further.

1/ Edit, if necessary, the template file catfile_template_cgbm.txt

   The setting PATH needs to be set to the top level location of the CGBM archive.

2/ Run the shell script "run_catfile_generation_cgbm.csh". This is a C-shell script to
   run the generic catalog file generator using the template for what files are in the
   archive (see above) and for the settings for the .cat file primary and 1st extenions.
   It needs to be run in the same directory where the templates are located and where the Perl script
   "make_generic_catfile.pl" is located.

   make_generic_catfile.pl is a Perl script to make "any mission" archive catalog files.
   Note that the logic for CGBM uses features introduce with v0.8 (8 Aug 2025) of this
   script: it is likely to not work with earlier versions.


## Contents

In this repository, you should find the following files:

`run_catfile_generation.csh`
A shell script example of running BurstCube .cat file generation for each of the 
* Transient/trigger data
* Daily data

`make_generic_catfile.pl`
The generic mission-independent code for generating catalogs

The following files are all templates used by make_generic_catfile.pl to make the four data types (above) for BurstCube. 

- cdfile_cat_basic.txt
- template_1stext_cat.txt
- template_primary_cat.txt
- catfile_template_bc_daily.txt
- catfile_template_bc_transient.txt


## Author

Jesse Allen (NASA/GSFC/HEASARC)

