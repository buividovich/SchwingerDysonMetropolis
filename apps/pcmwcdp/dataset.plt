reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace\\'
set pointsize 1.0
set bar 1.0
set style line 1  lt 1 lc rgb '#000000' lw 8 pt 5
set style line 2  lt 1 lc rgb '#FF0000' lw 8 pt 7
set style line 3  lt 1 lc rgb '#00FF00' lw 8 pt 9
set style line 4  lt 1 lc rgb '#0000FF' lw 8 pt 11
set style line 5  lt 1 lc rgb '#FF00FF' lw 8 pt 13
set style line 6  lt 1 lc rgb '#00FFFF' lw 8 pt 15
set style line 7  lt 1 lc rgb '#888888' lw 8 pt 1
set style line 11 lt 2 lc rgb '#000000' lw 8 pt 5
set style line 12 lt 2 lc rgb '#FF0000' lw 8 pt 7
set style line 13 lt 2 lc rgb '#00FF00' lw 8 pt 9
set style line 14 lt 2 lc rgb '#0000FF' lw 8 pt 11
set style line 15 lt 2 lc rgb '#FF00FF' lw 8 pt 13
set style line 16 lt 2 lc rgb '#00FFFF' lw 8 pt 15
set style line 17 lt 2 lc rgb '#888888' lw 8 pt 1

blueredcolortable = "#0000ff #1700e8 #2e00d1 #4600b9 #5d00a2 #74008b #8b0074 #a2005d #b90046 #d1002e #e80017 #ff0000"

plot_dir = 'G:\\LAT\\sd_metropolis\\plots\\pcm_wc_mspace'

set xlabel "1/m_{max}"
set yrange [-0.01:]
set key right bottom

set out sprintf("%s\\Gx_convergence_%s.eps", plot_dir, suffix)
set ylabel "G_x"
plot \
'Gx_'.suffix.'.mean' using ($1):($2):($3) notitle with yerrorbars ls 2

set out sprintf("%s\\mpnt_convergence_%s.eps", plot_dir, suffix)
set ylabel "G_{xy}(L_s/2)"
plot \
'mpnt_'.suffix.'.mean' using ($1):($2):($3) notitle with yerrorbars ls 2

set autoscale y
set out sprintf("%s\\link_convergence_%s.eps", plot_dir, suffix)
set ylabel "Link"
plot \
'link_'.suffix.'.mean' using ($1):($2):($3) notitle with yerrorbars ls 2

set out sprintf("%s\\correlator_%s.eps", plot_dir, suffix)
set yrange [-0.01:1.01]
set xrange [-0.2:54.2]
set xlabel "x"
set ylabel "G_{xy}(x)"
plot \
for [io=0:11] \
'Gxy_'.suffix.'.dat' using ($1):(column(2*io+2)):(column(2*io+3)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7

set autoscale y
set logscale y
set autoscale x
set key left bottom
set out sprintf("%s\\signs_%s.eps", plot_dir, suffix)
set ylabel "Mean sign"
set xlabel "Expansion order"
plot \
'signs_'.suffix.'.mean' using ($1):(abs($2)):($3) title 'Gx' with yerrorbars   ls 2,\
'signs_'.suffix.'.mean' using ($1):(abs($4)):($5) title 'Link' with yerrorbars ls 3
 

quit


