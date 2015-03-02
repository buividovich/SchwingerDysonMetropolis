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

OBS_FILES = system("ls G*.dat");
NFILES        = words(OBS_FILES)
FILE(i)       = word(OBS_FILES, i);
LABEL(i)      = FILE(i)[3:strstrt(FILE(i), ".dat")-1]

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\mean_G2_sign.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "Mean sign of G2"
plot \
for [i=1:NFILES] \
FILE(i) using ($1):($5) title LABEL(i) with points ls i 

unset logscale y
set yrange [0.0:0.05]

set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\observable_comparison_G20.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "G2(0)"
plot \
1/x notitle with lines ls 17, \
for [i=1:NFILES] \
FILE(i) using ($1):(abs($3-0.5)):($4) title LABEL(i) with yerrorbars ls i 

set out    'G:\\LAT\\sd_metropolis\\plots\\ssm\\observable_comparison_G21.eps'
set key right top
set xlabel "{/Symbol l}"
set ylabel "G2(1)"
plot \
1/x notitle with lines ls 17, \
for [i=1:NFILES] \
FILE(i) using ($1):(abs($7-0.5)):($8) title LABEL(i) with yerrorbars ls i 

#1/x notitle with lines ls 17, \
#1/x notitle with lines ls 17, \
