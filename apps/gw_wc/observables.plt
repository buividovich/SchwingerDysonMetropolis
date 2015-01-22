reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\gw_wc\\'
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

lc = 0.083333

set out    'G:\\LAT\\sd_metropolis\\plots\\gw_wc\\mean_recursion_depth.eps'
set logscale y
set key right top
set xlabel "{/Symbol l}"
set ylabel "Mean recursion depth"
plot\
'mc_stat_nmc5000000.dat'  using ($1):($5) title 'N_{mc} =  5x10^6'  with lines ls 2 

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\gw_wc\\mean_sign.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "Mean config. sign"
plot\
'mc_stat_nmc5000000.dat'  using ($1):($9) title 'N_{mc} =  5x10^6'  with lines ls 2 

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\gw_wc\\G2.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "G_n"
plot\
'G_nmc5000000.dat' using ($1):(1.0 - $2):($3)   title 'G_2'  with yerrorbars ls 2,\
0.25*x notitle with lines ls 3

