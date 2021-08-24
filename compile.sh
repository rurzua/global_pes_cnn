#!/bin/bash
#compile file work in macOS Catalina with MacPorts 2.7.1
gfortran-mp-9 -O3 -fdec-math -L/usr/lib -o globalpes.x src/myRKHS.f90 src/globalpes.f90
rm *.mod
