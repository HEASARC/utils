#!/bin/tcsh

./make_generic_catfile_v05.pl --cdfile cdfile_basic.txt \
  --primary template_primary_maxi.txt --extension template_1stext_maxi.txt \
  --range MJD --start 58800 --stop 59052 \
  --template catfile_template_maxi_v23.txt --compress --archive \
  |& tee maxi_catfile_58800_59052.log
