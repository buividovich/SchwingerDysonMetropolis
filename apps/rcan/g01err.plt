reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\rcan\\'
set pointsize 1.5
set bar 1
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

set out    'G:\\LAT\\sd_metropolis\\plots\\rcan\\G01err_gsc_l2.31_mplus.eps'
set key left bottom
set xlabel "m_0"
set ylabel "{/Symbol D}G_{01}/G_{01}"
set yrange [-0.5:0.5]
plot \
'mscan_l2.3100.dat' using ($1):($14) title 'm_{max} = 12' with lines ls 2,\
'mscan_l2.3100.dat' using ($1):($12) title 'm_{max} = 10' with lines ls 3,\
'mscan_l2.3100.dat' using ($1):($10) title 'm_{max} =  8' with lines ls 4,\
'mscan_l2.3100.dat' using ($1):($8)  title 'm_{max} =  6' with lines ls 5,\
0 notitle with lines ls 0 lw 1

quit

set out    'G:\\LAT\\sd_metropolis\\plots\\rcan\\G01err_gsc_l2.31_mplus.eps'
set key center top
set logscale y
#set yrange [-1.0:1.0]
set xlabel "1/(m_{max}+1)"
set ylabel "{/Symbol D}G_{01}/G_{01}"
plot \
'gsc_s2_l2.31_m0.00.dat' using ($2):(abs($5)) title 'm_0 = 0.0' with linespoints ls 2,\
'gsc_s2_l2.31_m0.10.dat' using ($2):(abs($5)) title 'm_0 = 0.1' with linespoints ls 3,\
'gsc_s2_l2.31_m0.20.dat' using ($2):(abs($5)) title 'm_0 = 0.2' with linespoints ls 4,\
'gsc_s2_l2.31_m0.50.dat' using ($2):(abs($5)) title 'm_0 = 0.5' with linespoints ls 5,\
'gsc_s2_l2.31_m1.00.dat' using ($2):(abs($5)) title 'm_0 = 1.0' with linespoints ls 6,\
'gsc_s2_l2.31_m2.00.dat' using ($2):(abs($5)) title 'm_0 = 2.0' with linespoints ls 7

fabs(x) = x

unset logscale y
set autoscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\rcan\\G01err_gsc_l2.31_mminus.eps'
set key right top
#set logscale y
set yrange [-0.1:0.1]
set xrange [0.0:0.15]
set xlabel "1/(m_{max}+1)"
set ylabel "{/Symbol D}G_{01}/G_{01}"
plot \
'gsc_s2_l2.31_m-0.6420.dat' using ($2):(fabs($5)) title 'm_0 = -0.642' with linespoints ls 1,\
'gsc_s2_l2.31_m-0.6440.dat' using ($2):(fabs($5)) title 'm_0 = -0.644' with linespoints ls 2,\
'gsc_s2_l2.31_m-0.6460.dat' using ($2):(fabs($5)) title 'm_0 = -0.646' with linespoints ls 3,\
'gsc_s2_l2.31_m-0.6480.dat' using ($2):(fabs($5)) title 'm_0 = -0.648' with linespoints ls 4,\
'gsc_s2_l2.31_m-0.65.dat' using ($2):(fabs($5)) title 'm_0 = -0.650' with linespoints ls 5,\
'gsc_s2_l2.31_m-0.6520.dat' using ($2):(fabs($5)) title 'm_0 = -0.652' with linespoints ls 6,\
'gsc_s2_l2.31_m-0.6540.dat' using ($2):(fabs($5)) title 'm_0 = -0.654' with linespoints ls 7,\
'gsc_s2_l2.31_m-0.6560.dat' using ($2):(fabs($5)) title 'm_0 = -0.656' with linespoints ls 8,\
'gsc_s2_l2.31_m-0.6580.dat' using ($2):(fabs($5)) title 'm_0 = -0.658' with linespoints ls 9,\
0 notitle with lines ls 1 lw 1

#'gsc_s2_l2.31_m-0.20.dat' using ($2):(fabs($5)) title 'm_0 = -0.20' with linespoints ls 2,\
#'gsc_s2_l2.31_m-0.25.dat' using ($2):(fabs($5)) title 'm_0 = -0.25' with linespoints ls 3,\
#'gsc_s2_l2.31_m-0.30.dat' using ($2):(fabs($5)) title 'm_0 = -0.30' with linespoints ls 4,\
#'gsc_s2_l2.31_m-0.35.dat' using ($2):(fabs($5)) title 'm_0 = -0.35' with linespoints ls 5,\
#'gsc_s2_l2.31_m-0.40.dat' using ($2):(fabs($5)) title 'm_0 = -0.40' with linespoints ls 6,\
#'gsc_s2_l2.31_m-0.45.dat' using ($2):(fabs($5)) title 'm_0 = -0.45' with linespoints ls 7,\
#'gsc_s2_l2.31_m-0.50.dat' using ($2):(fabs($5)) title 'm_0 = -0.50' with linespoints ls 8,\
