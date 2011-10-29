set  term X11
set auto
	filename(n) = sprintf("Statistical_data_concentration_X_%d", n)
	plot for [i=0:1] filename(i*1000).".sta" t "Number ".i w l
	



