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

a2(l) = (1 - sqrt(1 - 12*l))/(6*l)
G(n, l) = gamma(2*n + 1)/(gamma(n+1)*gamma(n+3))*a2(l)**n*(2*n + 2 - n*a2(l))

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\mean_G2_sign.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "Mean sign of G2"
plot\
'G_nmc20000000.dat'  using ($1):($6)  title 'G_2'  with lines ls 2,\
'G_nmc20000000.dat'  using ($1):($12) title 'G_4'  with lines ls 3

set logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\observable_comparison.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "G_n"
plot\
'G_nmc20000000.dat'  using ($1):($4):($5)    title 'G_2'  with yerrorbars ls 2,\
'G_nmc20000000.dat'  using ($1):($10):($11)  title 'G_4'  with yerrorbars ls 3,\
G(1, x)                                    notitle        with lines      ls 2,\
G(2, x)                                    notitle        with lines      ls 3,


quit

set logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\stack_stat.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "G_n"
plot \
for[i=1:8] \
'stack_stat_l'.sprintf("%2.4f", 0.01*i).'_nmc5000000.dat' using ($1):($2):($3)   title '{/Symbol l}='.sprintf("%2.4f", 0.01*i)  with yerrorbars ls i

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\ns_histories.eps'
set key left top
set xlabel "t_{MC}/1000"
set ylabel "m"
set yrange [0:60]
plot \
'ns_history_l0.1000_nmc10000.dat' using ($1/1000):($2+1) title '{/Symbol l}/{/Symbol l}_c=1.20' with lines ls 1 lw 2,\
'ns_history_l0.0850_nmc10000.dat' using ($1/1000):($2+1) title '{/Symbol l}/{/Symbol l}_c=1.02' with lines ls 2 lw 2,\
'ns_history_l0.0800_nmc10000.dat' using ($1/1000):($2+1) title '{/Symbol l}/{/Symbol l}_c=0.96' with lines ls 3 lw 2,\
'ns_history_l0.0500_nmc10000.dat' using ($1/1000):($2+1) title '{/Symbol l}/{/Symbol l}_c=0.60' with lines ls 4 lw 2


