reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\ssm\\'
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

PATSTR        = "G_smod";
OBS_FILES     = system("ls ".PATSTR."*.dat");
NFILES        = words(OBS_FILES)
FILE(i)       = word(OBS_FILES, i);
LABEL(i)      = FILE(i)[(strlen(PATSTR)+1):strstrt(FILE(i), ".dat")-1]

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\G0_smod.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "G(0)"
plot \
for [i=1:NFILES] \
FILE(i) using ($1):($4):($5) title "G(0), ".LABEL(i) with yerrorbars ls i,\
1 notitle with lines ls 0

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\G1_smod.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "G(1)"
plot \
for [i=1:NFILES] \
FILE(i) using ($1):($2):($3) title "G(1), ".LABEL(i) with yerrorbars ls i,\
1/x          notitle with lines ls 0 lc rgb '#000000',\
1/x + 1/x**3 notitle with lines ls 0 lc rgb '#FF0000'

#set logscale y
#set xrange [0:30]
#set yrange [0:15]
#set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\seq_len_hist.eps'
#set key right top
#set xlabel "n"
#set ylabel "c_n/c_{n+1}"
#plot \
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l80.0000.dat'  using 1:2:3 title "{/Symbol l} = 80"  with yerrorbars ls 1,\
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l85.0000.dat'  using 1:2:3 title "{/Symbol l} = 85"  with yerrorbars ls 2,\
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l90.0000.dat'  using 1:2:3 title "{/Symbol l} = 90"  with yerrorbars ls 3,\
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l95.0000.dat'  using 1:2:3 title "{/Symbol l} = 95"  with yerrorbars ls 4,\
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l100.0000.dat'  using 1:2:3 title "{/Symbol l} = 100"  with yerrorbars ls 5,\
#'stack_stat_smod_d2_s16_m-4.00_nmc50000000_l200.0000.dat'  using 1:2:3 title "{/Symbol l} = 200"  with yerrorbars ls 6
