#!/bin/bash

local_dir=/scr/animals1/Dropbox/01_code/topography/sandbox/individual_distance/0_behavdist/

while read i; do
  	echo "running ${i}"
  	mkdir /scr/animals1/lsd_freesurfer/${i}
  	mkdir /scr/animals1/lsd_freesurfer/${i}/surf
  	mkdir /scr/animals1/lsd_freesurfer/${i}/label
	cp /afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/${i}/surf/*.sphere.reg /scr/animals1/lsd_freesurfer/${i}/surf
	cp /afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/${i}/surf/*.smoothwm /scr/animals1/lsd_freesurfer/${i}/surf
	cp /afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/${i}/surf/*.pial /scr/animals1/lsd_freesurfer/${i}/surf
	cp /afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/${i}/label/*.aparc.a2009s.annot /scr/animals1/lsd_freesurfer/${i}/label
	cp /afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/${i}/label/*.cortex.label /scr/animals1/lsd_freesurfer/${i}/label

done<${1}


