reset
#set term postscript color enhanced landscape "Helvetica" 24
set term wxt
cd 'G:\\LAT\\sd_metropolis\\data\\hmm\\'
set pointsize 0.5
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

set logscale y
#set out    'G:\\LAT\\sd_metropolis\\plots\\hmm\\genus_test.eps'
set key right top
set xlabel "genus"
set ylabel "G_n^{(g)}"
plot\
'G_nmc20000000.dat'  using ($1):($2)   notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($4)   notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($6)   notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($8)   notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($10)  notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($12)  notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($14)  notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($16)  notitle   with points ls 2,\
'G_nmc20000000.dat'  using ($1):($3):(0.3)   notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($5):(0.3)   notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($7):(0.3)   notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($9):(0.3)   notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($11):(0.3)  notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($13):(0.3)  notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($15):(0.3)  notitle   with xerrorbars ls 4 lw 1,\
'G_nmc20000000.dat'  using ($1):($17):(0.3)  notitle   with xerrorbars ls 4 lw 1

