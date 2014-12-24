reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\HMM\\data\\'
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

IS(a, k) = (sqrt(pi)/2)*(a**(2+2*k))*gamma(k+1/2)/gamma(k+2)
a2(l)    = (1 - sqrt(1 - 12*l))/(6*l)
xm(l)    = 2*sqrt(a2(l))
G(l, n)  = IS(xm(l), n)*(1/2 - l*a2(l))/pi - l*IS(xm(l),n+1)/(2*pi)

lc = 0.083333

set out    'G:\\LAT\\sd_metropolis\\HMM\\mean_recursion_depth.eps'
set logscale y
set logscale x
set key right top
set xlabel "1 - {/Symbol l}/{/Symbol l}_c"
set ylabel "Mean recursion depth"
#set arrow from first 0.0,graph 0.0 to first 0.0, graph 1.0 nohead ls 7
plot\
'mc_stat_nmc200000.dat'   using ((lc-$1)/lc):($5) title 'N_{mc} =  2x10^5'  with lines ls 1 ,\
'mc_stat_nmc500000.dat'   using ((lc-$1)/lc):($5) title 'N_{mc} =  5x10^5'  with lines ls 2 ,\
'mc_stat_nmc1000000.dat'  using ((lc-$1)/lc):($5) title 'N_{mc} =  1x10^6'  with lines ls 3 ,\
'mc_stat_nmc2000000.dat'  using ((lc-$1)/lc):($5) title 'N_{mc} =  2x10^6'  with lines ls 4 ,\
'mc_stat_nmc5000000.dat'  using ((lc-$1)/lc):($5) title 'N_{mc} =  5x10^6'  with lines ls 5 

set out    'G:\\LAT\\sd_metropolis\\HMM\\mean_stack_size.eps'
set key right top
set xlabel "1 - {/Symbol l}/{/Symbol l}_c"
set ylabel "Mean stack size"
#set arrow from first 0.0,graph 0.0 to first 0.0, graph 1.0 nohead ls 7
plot\
 'mc_stat_nmc200000.dat'   using ((lc-$1)/lc):($6) title 'N_{mc} =  2x10^5'  with lines ls 1 ,\
 'mc_stat_nmc500000.dat'   using ((lc-$1)/lc):($6) title 'N_{mc} =  5x10^5'  with lines ls 2 ,\
'mc_stat_nmc1000000.dat'  using ((lc-$1)/lc):($6) title 'N_{mc} =  1x10^6'  with lines ls 3 ,\
'mc_stat_nmc2000000.dat'  using ((lc-$1)/lc):($6) title 'N_{mc} =  2x10^6'  with lines ls 4 ,\
'mc_stat_nmc5000000.dat'  using ((lc-$1)/lc):($6) title 'N_{mc} =  5x10^6'  with lines ls 5 

#unset logscale y
set out    'G:\\LAT\\sd_metropolis\\HMM\\ns_hist.eps'
set key right top
set xlabel "ns"
set ylabel "p(ns)"
plot\
'ns_hist_l0.0840_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0840'  with lines ls 1,\
'ns_hist_l0.0820_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0820'  with lines ls 2,\
'ns_hist_l0.0800_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0800'  with lines ls 3,\
'ns_hist_l0.0780_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0780'  with lines ls 4,\
'ns_hist_l0.0700_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0700'  with lines ls 5

#unset logscale y
set out    'G:\\LAT\\sd_metropolis\\HMM\\G_hist.eps'
set key right top
set xlabel "n"
set ylabel "G_n"
plot\
'G_l0.0840_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0840'  with lines ls 1,\
'G_l0.0820_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0820'  with lines ls 2,\
'G_l0.0800_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0800'  with lines ls 3,\
'G_l0.0780_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0780'  with lines ls 4,\
'G_l0.0700_nmc2000000.dat' using ($1):($2) title '{/Symbol l} = 0.0700'  with lines ls 5

set out    'G:\\LAT\\sd_metropolis\\HMM\\observable_comparison.eps'
set key right top
set xlabel "n"
set ylabel "G_n"
plot\
'G_l0.0820_nmc2000000.dat' using ($1):($2):($3) title '{/Symbol l} = 0.0820'  with yerrorbars ls 2

#/G(0.0820,$1)
unset logscale x
unset logscale y
set out    'G:\\LAT\\sd_metropolis\\HMM\\G_test.eps'
set key right top
set xlabel "{/Symbol l}"
set xrange [-lc:lc]
set ylabel "G_1({/Symbol l})"
plot\
 G(x,1) with lines ls 1

