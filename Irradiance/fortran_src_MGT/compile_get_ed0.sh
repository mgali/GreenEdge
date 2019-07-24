#!/bin/bash

# Remove old files
rm *.o
rm get_ed0

# Compile subroutines
gfortran -c interpol_ed0LUT_5nm_v2.f
gfortran -c read_ed0moins_LUT_5nm_v2.f 
gfortran -c read_ed0plus_LUT_5nm_v2.f 
gfortran -c sunpos.f

# Compile main code
gfortran get_ed0.f interpol_ed0LUT_5nm_v2.o read_ed0moins_LUT_5nm_v2.o read_ed0plus_LUT_5nm_v2.o sunpos.o -o get_ed0
