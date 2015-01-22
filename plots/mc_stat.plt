reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\'.app_name.'\\'
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

MC_STAT_FILES = system("ls mc_stat_nmc*.dat");
NFILES        = words(MC_STAT_FILES)
FILE(i)       = word(MC_STAT_FILES, i);
LABEL(i)      = substr(FILE(i), 11)

set out    'G:\\LAT\\sd_metropolis\\plots\\'.app_name.'\\mean_recursion_depth.eps'
set logscale y
set key right top
set xlabel "{/Symbol l}"
set ylabel "Mean recursion depth"
plot \
for [i=1:NFILES] \
FILE(i) using ($1):($5) with lines ls 2 

set out    'G:\\LAT\\sd_metropolis\\plots\\'.app_name.'\\mean_sign.eps'
set logscale y
set key right top
set xlabel "{/Symbol l}"
set ylabel "Mean sign (over all observables)"
plot \
for [i=1:NFILES] \
FILE(i) using ($1):($9) with lines ls 2 
