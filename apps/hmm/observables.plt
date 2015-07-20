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

set logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\G4.eps'
set key left top
set xlabel "{/Symbol l}"
set ylabel "G_4^{(g)}({/Symbol l})"
plot \
for [g=1:6] \
'G_nmc2000000000.dat'  using ($1):(column(3*g)):(column(3*g+1)) title 'g='.sprintf("%i", g) with yerrorbars ls g,\
for [g=1:6] \
'hmm_G4_analytics.dat' using ($1):(column(g+1)):(1.0)  every ::1          notitle         with lines  ls g

unset logscale y
set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\G4_ratios.eps'
set key left bottom
set xrange [0.0:1.0]
set xlabel "{/Symbol l}/{/Symbol l}_c"
set ylabel "G_4(MC)/G_4(analytics)"
plot \
1 with lines notitle ls 1, \
for [g=0:6] \
'rg'.sprintf("%i", g).'.dat' using (12*$1):($2):($3) title 'g='.sprintf("%i", g) with yerrorbars ls (g+1)


quit

