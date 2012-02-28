#!/bin/sh
gnuplot -persist<<EOF




set multiplot

set origin 0.0,0.0
set size 0.333,0.333
plot "$1Velocity_ave_max.txt" u 1 w l title "v_ave"



set origin 0.333,0.0
set size 0.333,0.333
plot "$1Velocity_ave_max.txt" u 2 w l title "v_max"


set origin 0.666,0.0
set size 0.333,0.333
plot "$1bodyforce.txt" u 1 w l title "saturation of 1"



set origin 0.0,0.333
set size 0.333,0.333
plot "$1Relative_Permeability_Component1.txt" u 5 w l title "pressure left",\
"$1Relative_Permeability_Component1.txt" u 6 w l title "pressure right"
#"$1Relative_Permeability_Component2.txt" u 3 title "LOCAL_Rel_Perm_Both"
#plot "$1bodyforce.txt" u 2 w l title "saturation of 2"


set origin 0.333,0.333
set size 0.333,0.333
plot "$1Relative_Permeability_Component1.txt" u 2:3 w l title "LOCAL_Rel_Perm1",\
"$1Relative_Permeability_Component2.txt" u 4:3 w l title "LOCAL_Rel_Perm2"



set origin 0.666,0.333
set size 0.333,0.333
plot "$1Relative_Permeability_Component1.txt" u 2:1 w l title "Rel_Perm1",\
"$1Relative_Permeability_Component2.txt" u 4:1 w l title "Rel_Perm2"


set origin 0.333,0.666
set size 0.333,0.333
plot "$1Relative_Permeability_Component1.txt" u 4:3 w l title "LOCAL_Rel_Perm1",\
"$1Relative_Permeability_Component2.txt" u 2:3 w l title "LOCAL_Rel_Perm2"



set origin 0.666,0.666
set size 0.333,0.333
plot "$1Relative_Permeability_Component1.txt" u 4:1 w l title "Rel_Perm1",\
"$1Relative_Permeability_Component2.txt" u 2:1 w l title "Rel_Perm2"

set origin 0.0,0.666
set size 0.333,0.333
plot "$1bodyforce.txt" u 1 w l title "saturation of 1",\
"$1bodyforce.txt" u 2 w l title "saturation of 2"

pause -1
