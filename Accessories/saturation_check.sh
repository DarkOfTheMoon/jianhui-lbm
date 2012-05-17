#!/bin/sh
gnuplot -persist<<EOF




set multiplot

set origin 0.0,0.0
set size 0.5,1.0
plot "$1bodyforce.txt" u 1 w l title "saturation of 1"


set origin 0.5,0.0
set size 0.5,1.0
set logscale y
plot "$1bodyforce.txt" u 3  title "saturation error"


pause -1
