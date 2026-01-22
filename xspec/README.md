# Xspec Utilities

Various utilities useful to users and maintainers of [Xspec](https://heasarc.gsfc.nasa.gov/docs/software/xspec/index.html).

## Contents

`make_spex_files_for_xspec.py` - Program to create apec-style files by running a spex model at various temperature grid points. The lines and continuum are then gathered up, and stored in an apec-format table model for xspec. These files are distributed with [spectral model data files](https://heasarc.gsfc.nasa.gov/docs/software/xspec/modeldata.html).

`xspec_utils.py` - a set of pyxpsec functions to read xcm files created by the command-line version of xspec, return lists of line identifications, and other tasks

