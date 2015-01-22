reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\hmm\\'
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

set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\mean_recursion_depth.eps'
set logscale y
set key right top
set xlabel "1 - |{/Symbol l}|/{/Symbol l}_c"
set ylabel "Mean recursion depth"
#set arrow from first 0.0,graph 0.0 to first 0.0, graph 1.0 nohead ls 7
plot\
'mc_stat_nmc5000000.dat'  using ((lc-abs($1))/lc):($5) title 'N_{mc} =  5x10^6'  with lines ls 2 

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\mean_stack_size.eps'
set key right top
set xlabel "1 - |{/Symbol l}|/{/Symbol l}_c"
set ylabel "Mean stack size"
#set arrow from first 0.0,graph 0.0 to first 0.0, graph 1.0 nohead ls 7
plot\
'mc_stat_nmc5000000.dat'  using ((lc-abs($1))/lc):($6) title 'N_{mc} =  5x10^6'  with lines ls 2 

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\mean_G2_sign.eps'
set key right top
set xlabel "1 - |{/Symbol l}|/{/Symbol l}_c"
set ylabel "Mean sign of G2"
#set arrow from first 0.0,graph 0.0 to first 0.0, graph 1.0 nohead ls 7
plot\
'G_nmc5000000.dat'  using ((lc-abs($1))/lc):($4) title 'N_{mc} =  5x10^6'  with lines ls 2 

set logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\observable_comparison.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "G_n"
plot\
'G_nmc5000000.dat' using ($1):($2):($3)   title 'G_2'  with yerrorbars ls 2,\
'G_nmc5000000.dat' using ($1):($5)        notitle      with lines      ls 3,\
'G_nmc5000000.dat' using ($1):($6):($7)   title 'G_4'  with yerrorbars ls 4,\
'G_nmc5000000.dat' using ($1):($9)        notitle      with lines      ls 5,\
'G_nmc5000000.dat' using ($1):($10):($11) title 'G_6'  with yerrorbars ls 6,\
'G_nmc5000000.dat' using ($1):($13)       notitle      with lines      ls 7

