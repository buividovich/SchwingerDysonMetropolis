reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\stereo\\'
set pointsize 1.5
set bar 1
set style line 1  lt 1 lc rgb '#000000' lw 8 pt 5
set style line 2  lt 2 lc rgb '#FF0000' lw 8 pt 7
set style line 3  lt 3 lc rgb '#00FF00' lw 8 pt 9
set style line 4  lt 4 lc rgb '#0000FF' lw 8 pt 11
set style line 5  lt 5 lc rgb '#FF00FF' lw 8 pt 13
set style line 6  lt 6 lc rgb '#00FFFF' lw 8 pt 15
set style line 7  lt 7 lc rgb '#888888' lw 8 pt 1
set style line 8  lt 8 lc rgb '#88FFFF' lw 8 pt 1
set style line 11 lt 2 lc rgb '#000000' lw 8 pt 5
set style line 12 lt 2 lc rgb '#FF0000' lw 8 pt 7
set style line 13 lt 2 lc rgb '#00FF00' lw 8 pt 9
set style line 14 lt 2 lc rgb '#0000FF' lw 8 pt 11
set style line 15 lt 2 lc rgb '#FF00FF' lw 8 pt 13
set style line 16 lt 2 lc rgb '#00FFFF' lw 8 pt 15
set style line 17 lt 2 lc rgb '#888888' lw 8 pt 1

set out    'G:\\LAT\\sd_metropolis\\doc\\plots\\Gx.eps'
set key left top
set xlabel "1/max.order"
set ylabel "<Gx>"
set yrange [0:]
set xrange [0:]
plot \
'Gx.dat' using (1/($1+1)):($2):($3) notitle with yerrorbar ls 2 lw 2

set out    'G:\\LAT\\sd_metropolis\\doc\\plots\\Gx_sign.eps'
set key left top
set xlabel "1/max.order"
set ylabel "<Gx>"
set autoscale y
set logscale y
set xrange [0:]
plot \
'Gx.dat' using (1/($1+1)):(abs($4)) notitle with points ls 2 lw 2

set out    'G:\\LAT\\sd_metropolis\\doc\\plots\\Gxy.eps'
set key left top
set xlabel "p"
set ylabel "G(p)"
set autoscale y
set logscale y
set xrange [0:]
set yrange [0.005:1.0]
plot \
'Gxy_o0.dat' using ($1):($2):($3) title 'order 0' with yerrorbars ls 2 lw 2,\
'Gxy_o1.dat' using ($1):($2):($3) title 'order 1' with yerrorbars ls 3 lw 2,\
'Gxy_o2.dat' using ($1):($2):($3) title 'order 2' with yerrorbars ls 4 lw 2,\
'Gxy_o3.dat' using ($1):($2):($3) title 'order 3' with yerrorbars ls 5 lw 2,\
'Gxy_o4.dat' using ($1):($2):($3) title 'order 4' with yerrorbars ls 6 lw 2,\
'Gxy_o5.dat' using ($1):($2):($3) title 'order 5' with yerrorbars ls 7 lw 2
