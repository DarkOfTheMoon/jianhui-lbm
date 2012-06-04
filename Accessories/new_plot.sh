#!/bin/sh
gnuplot -persist<<EOF
set multiplot
set origin 0.0,0.0
set size 0.5,0.5
plot "$1Relative_Permeability_Component1.txt" u 1 w l title "Rel_Perm1",\
"$1Relative_Permeability_Component2.txt" u 1 w l title "Rel_Perm2"

set origin 0.5,0.0
set size 0.5,0.5
plot "$1new.dat" u 1 w l title "Least_Rel_Perm1",\
"$1new.dat" u 2 w l title "Least_Rel_Perm2"

set origin 0.0,0.5
set size 0.5,0.5
plot "$1new.dat" u 3 w l title "AVE_Rel_Perm1",\
"$1new.dat" u 4 w l title "AVE_Rel_Perm2"

set origin 0.5,0.5
set size 0.5,0.5
set log y
plot "$1new.dat" u 5 w l title "DEV_Rel_Perm1",\
"$1new.dat" u 6 w l title "DEV_Rel_Perm2"
unset log y






