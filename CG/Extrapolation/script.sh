#!/bin/bash

cat raw.table | awk '{print $2, $3}' > u_bumper.dat
cat raw.table | awk '{print $2, $4}' > forces_bumper.dat

python extrapolate.py
python extrapolate_2.py

gfortran -O3 pot_table.f90 -o table.exe
./table.exe

rm *dat*
