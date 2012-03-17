#!/bin/sh
gnuplot -persist<<EOF




set multiplot

set origin 0.0,0.0
set size 0.5,1.0
#plot "$1prosity_output.dat" u 1:2  title "porosity of subdomains in different sizes"
plot "$1prosity_output.dat_sta" u 1:5  w l title "Sandard deviation of subdomians with different sizes"


set origin 0.5,0.0
set size 0.5,1.0

plot "$1prosity_output.dat_sta" u 1:2:3:4 with yerrorbars title "error bas of subdomains",\
"$1prosity_output.dat_sub1" title "subdomains with increasing sizes"

pause -1
