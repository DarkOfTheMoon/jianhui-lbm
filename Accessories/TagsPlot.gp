#set terminal png         # gnuplot recommends setting terminal before output
#set output "output.png"  # The output filename; to be set after setting

#set ylabel "$sin(\\theta)$"  #use two // to represent / in latex

set terminal png         # gnuplot recommends setting terminal before output
set output "output.png"  # The output filename; to be set after setting

#-0.1*x+0.1*y+0.505 "TAGSDATA.dat" u 13:15:14 t "tire location" w l,\
#set auto
#set xrange[-100:450]
#set yrange[-100:450]
#set zrange[0:30]

#set auto
#splot\
#"TAGSDATA.dat" u 1:3:2 t "wheelFL position" w l,\
#"TAGSDATA.dat" u 4:6:5 t "wheelFR position" w l,\
#"TAGSDATA.dat" u 10:12:11  t "wheelBL position" w l,\
#"TAGSDATA.dat" u 16:18:17  t "wheelBR position" w l,\
#"TAGSDATA.dat" u 7:9:8  t "wheelML position" w l,\
#"TAGSDATA.dat" u 13:15:14  t "wheelMR position" w l

#y(x) = a*x**2 + b*x + c
#fit y(x) "data.txt" via a,b,c
#plot "data.txt", y(x)

set xrange [-pi:pi]                       # we want only one cycle
set xtics ("0" 0, \
	        "90" pi/2, "-90" -pi/2, \
			"" pi/4 1, "" -pi/4 1,  \
			"" 3*pi/4 1, "" -3*pi/4 1)
set grid
set xlabel "Angle,\n in degrees"
set ylabel "sin(angle)"
plot sin(x)

#set term post enh  		 # enhanced PostScript, essentially PostScript
 		 		 # with bounding boxes
#set out 'gplt.eps'
#set xlabel '{/Symbol q_1}
#set ylabel 'sin^2({/Symbol q_1})'
#plot sin(x)**2

#plot "data.txt" using (log10($1)):(log10($2))

#y1(x) = A1*x**b1
#y2(x) = A2*x**b2
#fit [0.01:0.3] y1(x) "data.txt" via A1,b1
#fit [3:30] y2(x) "data.txt" via A2,b2
#plot "data.txt", y1(x), y2(x)

