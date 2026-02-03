#!/bin/tcsh
#
# C-shell script to generate .cat files for CALET CGBM data

./make_generic_catfile.pl -cdfile cdfile_cat_basic.txt \
  --primary template_primary_cat.txt --extension template_1stext_cat.txt \
<<<<<<< HEAD
  --template catfile_template_cgbm.txt --compress --archive 

=======
  --template catfile_template_cgbm.txt --compress --archive 
>>>>>>> 65cca133d348eee9c235b9f3b71b51b61db90ddf
