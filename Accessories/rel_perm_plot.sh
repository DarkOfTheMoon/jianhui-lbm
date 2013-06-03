#!/bin/sh
gnuplot -persist<<EOF


set multiplot


set origin 0.0,0.0
set size 0.5,1.0
plot [0:1][0:1] "berea_shell_dri_nw.dat" u 1:2 w l lc 1 lw 2  title "oil",\
"berea_shell_dri_w.dat" u 1:2 w l lc 3 lw 2  title "water",\
"Summary_Ave.outdat" u 1:2 w p lc 5 lw 2 pt 7 ps 3 title "LB water",\
"Summary_Ave.outdat" u 1:3 w p lc 7 lw 2 pt 5 ps 3 title "LB oil"

set origin 0.5,0.0
set size 0.5,1.0
plot [0:1][0:1] "berea_shell_dri_nw.dat" u 1:2 w l lc 1 lw 2  title "oil",\
"berea_shell_dri_w.dat" u 1:2 w l lc 3 lw 2  title "water",\
"Summary_Least_Square.outdat" u 1:2 w p lc 5 lw 2 pt 7 ps 3 title "LB water least",\
"Summary_Least_Square.outdat" u 1:3 w p lc 7 lw 2 pt 5 ps 3 title "LB oil least"


#plot [0:1][0:1] "berea_shell_dri_nw.dat" u 1:2 w p lc 1 lw 2 pt 9 ps 2 title "oil",\
#"berea_shell_dri_w.dat" u 1:2 w p lc 3 lw 2 pt 11 ps 2  title "water",\
#"Summary_Ave.outdat" u 1:2 w p lc 5 lw 2 pt 7 ps 3 title "LB water",\
#"Summary_Ave.outdat" u 1:3 w p lc 7 lw 2 pt 5 ps 3 title "LB oil"


#set origin 0.0,0.666
#set size 0.333,0.333
#plot "$1bodyforce.txt" u 1 w l title "saturation of 1",\
#"$1bodyforce.txt" u 2 w l title "saturation of 2"
#plot "$1Capillary_Pressure.txt" u 1:3 w lp title "Capillary pressure phase1"


pause -1
