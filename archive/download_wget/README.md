# Downloading Data

## Introduction

This script provides a convenient way to download a requested set of data from the HEASARC archives.

See the [HEASARC webpage about command line downloads](https://heasarc.gsfc.nasa.gov/docs/cookbook/command_line_downloads.html) for more information.


## Usage 

```
download_wget <url> 
```

Downloads data to the local computer from the HEASARC archive using wget command and a url. 

The possible type of downloads are:

1. single directory corresponding to a specific observation
2. range of directories corresponding to a range of observations 
3. single file 
4. type of files 

If the transfer is interruped, the user may restart the same command in the same directory and only
the remaining files are downloaded.  

### 1) single directory corresponding to a specific observation

The command for a single directory is : 

```
download_wget.pl https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/
```

On the local computer the downloaded data are in the directory 1050020180/ maintaining  
  with the archive structure.

### 2. range of directories corresponding to a range of observations 

The command for a range of directories is :

```
download_wget.pl "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/105002018[01]"
```

On the local computer the downloaded data are in the directories 1050020180/ and 
1050020181/. The archive structure is maintained within each directory.

### 3. single file 

The command for a single file is :

```
download_wget.pl https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/1050020180/auxil/ni1050020180.att.gz
```

On the local computer the file is downloaded in the directory where the script is invoked. 

### 4. type of files 

The command for multiple files with the same pattern : 

```
download_wget.pl "https://heasarc.gsfc.nasa.gov/FTP/nicer/data/obs/2018_01/*/auxil/ni*.att.gz"
```

On the local computer the files are downloaded in the directory where the script is invoked.

The archive structure is not maintained.

Within a given url wildcards are allowed to search specific files that are located in different 
directories with a well defined pattern. The wildcard completion allowed are * and [ ]. A 
maximum of two non consegutive wildcards are allowed for url.


