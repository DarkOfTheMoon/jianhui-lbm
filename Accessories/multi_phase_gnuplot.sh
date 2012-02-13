#!/bin/sh
gnuplot -persist<<EOF




set multiplot

set origin 0.0,0.0
set size 0.333,0.5
plot "$1Velocity_ave_max.txt" u 1 w l title "v_ave"



set origin 0.333,0.0
set size 0.333,0.5
plot "$1Velocity_ave_max.txt" u 2 w l title "v_max"


set origin 0.666,0.0
set size 0.333,0.5
plot "$1bodyforce.txt" u 1 w l title "saturation of 1"



set origin 0.0,0.5
set size 0.333,0.5
plot "$1Relative_Permeability_Component1.txt" u 1 w l title "Rel_Perm_1"



set origin 0.333,0.5
set size 0.333,0.5
plot "$1Relative_Permeability_Component1.txt" u 3 w l,\
"$1Relative_Permeability_Component2.txt" u 3 w l title "LOCAL_Rel_Perm_Both"



set origin 0.666,0.5
set size 0.333,0.5
plot "$1Relative_Permeability_Component1.txt" u 1 w l,\
"$1Relative_Permeability_Component2.txt" u 1 w l title "Rel_Perm_Both"
pause -1
