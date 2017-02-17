#bin/bash

cd lib
ifort -c qsort_c.f90 -O2
ifort -c pointers.f90 -O2
ifort -c sg_gridutils.f90 -O2
ifort -c sg_griders.f90 -O2
ifort -c sg_utils.f90 -O2
cd ..
ifort -o subgrider main/sg_main.f90 lib/qsort_c.f90 lib/pointers.f90 lib/sg_gridutils.f90 lib/sg_griders.f90 lib/sg_utils.f90 -O2 -module lib

