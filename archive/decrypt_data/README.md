# Decrypting Data

## Introduction

Select observational data in the archive are temporarily encrypted during a limited-use period. Each file can be decrypted one at a time. This script is a convenient way to loop through a directory and decrypt the entire contents in one command.

For more information, see the [HEASARC cookbook pages](https://heasarc.gsfc.nasa.gov/docs/cookbook/).

## Description

This script decrypts all *.pgp or *.gpg files in a directory tree
specified with the -d option.  Individual files may be decrypted
by giving them as arguments to the script.

The script tests if you have PGP and GPG installed.
By default the script uses PGP, however, for any of the following cases
GPG is used:

  * The `-g` option is set
  * The files end in `.gpg`
  * PGP is unavailable

The script has to work at least in these platforms and with these
versions of PGP and GPG: 

- For Linux and Sun, it works with PGP Version 6.5.8.(Command='pgp')
- For Linux, it works with GPG (GnuPG) Version 1.0.7 (Command='gpg')
- For Macintosh, it works with GPG (GnuPG) Version 1.2.4 (Command='gpg')
                                  (GnuPG/MacGPG2) Vaersion 2.0.22  

The extension of the files to be decrypted has to be `.pgp` or `.gpg`

Options:

```
-h            This help
-v            Version number
-p password   Set password used to decrypt
-d directory  Set directory to recursively decrypt
-r            Remove the encrypted files *.pgp or *.gpg
-f            Force the use of PGP
-g            Force the use of GPG
-q            Quietly decrypt
```

Suggested usage:

```
> decrypt_data.pl -d directory

Using GPG 1.2.6

Enter the password:
```

Then, enter the password or decrypting key at the prompt. 

