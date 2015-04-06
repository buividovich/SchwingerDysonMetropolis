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

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\acceptance_vs_pplus.eps'
set key right top
set xlabel "p_+"
set ylabel "Acceptance rate"
plot\
'mc_stat_acctest_nmc20000000_l0.0400.dat'  using ($12):($7) title 'N_{mc} =  2x10^7, {/Symbol l} = 0.04'  with lines ls 2,\
'mc_stat_acctest_nmc20000000_l0.0600.dat'  using ($12):($7) title 'N_{mc} =  2x10^7, {/Symbol l} = 0.06'  with lines ls 3,\
'mc_stat_acctest_nmc20000000_l0.0800.dat'  using ($12):($7) title 'N_{mc} =  2x10^7, {/Symbol l} = 0.08'  with lines ls 4
