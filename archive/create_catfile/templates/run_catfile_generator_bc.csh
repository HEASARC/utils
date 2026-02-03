#!/bin/tcsh
#
# C-shell script to generate .cat files for BurstCube daily and trigger data

make_generic_catfile.pl -cdfile cdfile_cat_basic.txt \
  --primary template_primary_cat.txt --extension template_1stext_cat.txt \
  --template catfile_template_bc_daily.txt --compress --archive 

make_generic_catfile.pl -cdfile cdfile_cat_basic.txt \
  --primary template_primary_cat.txt --extension template_1stext_cat.txt \
  --template catfile_template_bc_transient.txt --compress --archive

