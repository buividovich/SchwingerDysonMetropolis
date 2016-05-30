reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\smlm\\'
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

fy(y)    = y
dy(y, e) = e

set out    'G:\\LAT\\sd_metropolis\\plots\\smlm\\G2_total_vs_ao_l1.20.eps'
set key center top
set xlabel "ao"
set yrange [-0.8:0.8]
set xrange [0.0:6.5]
set ylabel "G2_total"
plot \
0 with lines ls 1,\
'G2_nmc30000000_a0.29_l1.20_s2.dat'    using ($1+0.0):(fy($2)):(dy($2,$3)) title 'L_s = 2'    with yerrorbars ls 2,\
'G2_nmc30000000_a0.29_l1.20_s20.dat'   using ($1+0.1):(fy($2)):(dy($2,$3)) title 'L_s = 20'   with yerrorbars ls 3,\
'G2_nmc30000000_a0.29_l1.20_s200.dat'  using ($1+0.2):(fy($2)):(dy($2,$3)) title 'L_s = 200'  with yerrorbars ls 4,\
'G2_nmc30000000_a0.29_l1.20_s2000.dat' using ($1+0.3):(fy($2)):(dy($2,$3)) title 'L_s = 2000' with yerrorbars ls 5

set out    'G:\\LAT\\sd_metropolis\\plots\\smlm\\G2_total_vs_ao_l2.30.eps'
set key center top
set xlabel "ao"
set yrange [-0.8:0.8]
set xrange [0.0:6.5]
set ylabel "G2_total"
plot \
0 with lines ls 1,\
'G2_nmc30000000_a0.25_l2.30_s2.dat'    using ($1+0.0):(fy($2)):(dy($2,$3)) title 'L_s = 2'    with yerrorbars ls 2,\
'G2_nmc30000000_a0.25_l2.30_s20.dat'   using ($1+0.1):(fy($2)):(dy($2,$3)) title 'L_s = 20'   with yerrorbars ls 3,\
'G2_nmc30000000_a0.25_l2.30_s200.dat'  using ($1+0.2):(fy($2)):(dy($2,$3)) title 'L_s = 200'  with yerrorbars ls 4,\
'G2_nmc30000000_a0.25_l2.30_s2000.dat' using ($1+0.3):(fy($2)):(dy($2,$3)) title 'L_s = 2000' with yerrorbars ls 5

set out    'G:\\LAT\\sd_metropolis\\plots\\smlm\\G2_total_vs_ao_l0.60.eps'
set key center top
set xlabel "ao"
set yrange [-0.8:0.8]
set xrange [0.0:6.5]
set ylabel "G2_total"
plot \
0 with lines ls 1,\
'G2_nmc30000000_a0.20_l0.60_s2.dat'    using ($1+0.0):(fy($2)):(dy($2,$3)) title 'L_s = 2'    with yerrorbars ls 2,\
'G2_nmc30000000_a0.20_l0.60_s20.dat'   using ($1+0.1):(fy($2)):(dy($2,$3)) title 'L_s = 20'   with yerrorbars ls 3,\
'G2_nmc30000000_a0.20_l0.60_s200.dat'  using ($1+0.2):(fy($2)):(dy($2,$3)) title 'L_s = 200'  with yerrorbars ls 4,\
'G2_nmc30000000_a0.20_l0.60_s2000.dat' using ($1+0.3):(fy($2)):(dy($2,$3)) title 'L_s = 2000' with yerrorbars ls 5

set out    'G:\\LAT\\sd_metropolis\\plots\\smlm\\G2_total_vs_ao_l0.60_morestat.eps'
set key center top
set xlabel "ao"
set yrange [-1.0:1.0]
set xrange [0.0:9.5]
set ylabel "G2_total"
plot \
0 with lines ls 1,\
'G2_nmc1500000000_a0.20_l0.60_s2.dat'    using ($1+0.0):(fy($2)):(dy($2,$3)) title 'L_s = 2'    with yerrorbars ls 2,\
'G2_nmc1500000000_a0.20_l0.60_s20.dat'   using ($1+0.1):(fy($2)):(dy($2,$3)) title 'L_s = 20'   with yerrorbars ls 3,\
'G2_nmc1500000000_a0.20_l0.60_s200.dat'  using ($1+0.2):(fy($2)):(dy($2,$3)) title 'L_s = 200'  with yerrorbars ls 4,\
'G2_nmc1500000000_a0.20_l0.60_s2000.dat' using ($1+0.3):(fy($2)):(dy($2,$3)) title 'L_s = 2000' with yerrorbars ls 5
