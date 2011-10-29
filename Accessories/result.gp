#set  term X11
#set auto
#plot\
#"DoublePendulumResult.dat" u 8 t "x" w l,\
#"DoublePendulumResult.dat" u 9 t "y" w l,\
#"DoublePendulumResult.dat" u 10 t "z" w l
#plot\
#"DoublePendulumResult.dat" u 11 t "pivot FL" w l,\
#"DoublePendulumResult.dat" u 12 t "PIVOT FR" w l,\
#"DoublePendulumResult.dat" u 16 t "PIVOT ML" w l,\
#"DoublePendulumResult.dat" u 17 t "pivot BL" w l,\
#"DoublePendulumResult.dat" u 21 t "PIVOT MR" w l,\
#"DoublePendulumResult.dat" u 22 t "PIVOT BR" w l

#---------position--------------
#set size 0.33, 0.5
#set origin 0.33, 0.5
#plot\
#"DoublePendulumResult.dat" u 2 t "position x" w l,\
#"DoublePendulumResult.dat" u 3 t "position y" w l,\
#"DoublePendulumResult.dat" u 4 t "position z" w l


set size 1,1
set origin 0,0
set multiplot
set auto

set size 0.25, 0.33
set origin 0.0, 0.66
set yrange[-5:5]
plot\
"DoublePendulumResult.dat" u 26 t "yConstrain 0" w l,\
"DoublePendulumResult.dat" u 27 t "yConstrain 1" w l,\
"DoublePendulumResult.dat" u 28 t "yConstrain 2" w l,\
"DoublePendulumResult.dat" u 29 t "yConstrain 3" w l,\
"DoublePendulumResult.dat" u 30 t "yConstrain 4" w l,\
"DoublePendulumResult.dat" u 31 t "yConstrain 5" w l

set auto
set size 0.25, 0.33
set origin 0.25, 0.66
plot\
"DoublePendulumResult.dat" u 2 t "position x" w l,\
"DoublePendulumResult.dat" u 3 t "position y" w l,\
"DoublePendulumResult.dat" u 4 t "position z" w l


set auto
set size 0.25, 0.33
set origin 0.5, 0.66
plot\
"DoublePendulumResult.dat" u 5 t "x rotation degree" w l,\
"DoublePendulumResult.dat" u 6 t "y rotation degree" w l,\
"DoublePendulumResult.dat" u 7 t "z rotation degree" w l

set size 0.25, 0.33
set origin 0.75, 0.66

plot\
"DoublePendulumResult.dat" u 8 t "Axil F" w l,\
"DoublePendulumResult.dat" u 13 t "Axil BL" w l,\
"DoublePendulumResult.dat" u 18 t "Axil BR" w l

set auto
set size 0.25, 0.33
set origin 0.0, 0.33
plot\
"DoublePendulumResult.dat" u 9 t "pivot FL (PID)" w l,\
"DoublePendulumResult.dat" u 10 t "PIVOT FR (PID)" w l,\
"DoublePendulumResult.dat" u 14 t "PIVOT ML (PID)" w l,\
"DoublePendulumResult.dat" u 15 t "pivot BL (PID)" w l,\
"DoublePendulumResult.dat" u 19 t "PIVOT MR (PID)" w l,\
"DoublePendulumResult.dat" u 20 t "PIVOT BR (PID)" w l


set size 0.25,0.33
set origin 0.25,0.33
plot\
"DoublePendulumResult.dat" u 23 t "velocity x" w l,\
"DoublePendulumResult.dat" u 24 t "velocity y" w l,\
"DoublePendulumResult.dat" u 25 t "velocity z" w l,\
"DoublePendulumResult.dat" u ($23*$23+$24*$24+$25*$25)**0.5 t "velocity average" w l


set size 0.25, 0.33
set origin 0.5,0.33
plot\
"DoublePendulumResult.dat" u 11 t "wheel FL" w l,\
"DoublePendulumResult.dat" u 12 t "wheel FR" w l,\
"DoublePendulumResult.dat" u 16 t "wheel ML" w l,\
"DoublePendulumResult.dat" u 17 t "wheel BL" w l,\
"DoublePendulumResult.dat" u 21 t "wheel MR" w l,\
"DoublePendulumResult.dat" u 22 t "wheel BR" w l

set yrange[-5:300]
set size 0.25,0.33
set origin 0.75,0.33
plot\
"DoublePendulumResult.dat" u 32 t "ReactionForce 0" w l

set auto
set size 0.25,0.33
set origin 0.0,0.0
plot\
"DoublePendulumResult.dat" u 33 t "Wheel R Velocity FL" w l,\
"DoublePendulumResult.dat" u 34 t "Wheel R Velocity FR" w l,\
"DoublePendulumResult.dat" u 35 t "Wheel R Velocity ML" w l,\
"DoublePendulumResult.dat" u 36 t "Wheel R Velocity BL" w l,\
"DoublePendulumResult.dat" u 37 t "Wheel R Velocity MR" w l,\
"DoublePendulumResult.dat" u 38 t "Wheel R Velocity BR" w l


set size 0.25,0.33
set origin 0.25,0.0
set auto
splot\
"TAGSDATA.dat" u 1:3:2 t "wheelFL position" w l,\
"TAGSDATA.dat" u 4:6:5 t "wheelFR position" w l,\
"TAGSDATA.dat" u 10:12:11  t "wheelBL position" w l,\
"TAGSDATA.dat" u 16:18:17  t "wheelBR position" w l,\
"TAGSDATA.dat" u 7:9:8  t "wheelML position" w l,\
"TAGSDATA.dat" u 13:15:14  t "wheelMR position" w l

set size 0.25, 0.33
set origin 0.50, 0.0
set yrange[-5:5]
plot\
"DoublePendulumResult.dat" u 39 t "y2Constrain 0" w l,\
"DoublePendulumResult.dat" u 40 t "y2Constrain 1" w l,\
"DoublePendulumResult.dat" u 41 t "y2Constrain 2" w l,\
"DoublePendulumResult.dat" u 42 t "y2Constrain 3" w l,\
"DoublePendulumResult.dat" u 43 t "y2Constrain 4" w l,\
"DoublePendulumResult.dat" u 44 t "y2Constrain 5" w l


  
unset multiplot


