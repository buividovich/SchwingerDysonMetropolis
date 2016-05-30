reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\rcan\\'
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

n = 3.0
fx(x) = x #(x<0? -1 : +1)*(abs(x)**(1/n))
fy(y) = y #(y<0? -1 : +1)*(abs(y)**(1/n))

set out    'G:\\LAT\\sd_metropolis\\plots\\rcan\\stereo_G01_conv.eps'
set key left top
set xlabel "1/m_{max}"
set ylabel "{/Symbol D}G_{01}/G_{01}"
set yrange [-0.1:0.1]
set xrange [0.0:1.0]
plot \
0 notitle with lines ls 0 lw 1,\
'stereo_conv_l2.3100.dat'  using (fx($2)):(fy($4)) title '{/Symbol l}=0.50'  with linespoints ls 2,\
'stereo_conv_l2.3100.dat'  using (fx($2)):(fy($4)) title '{/Symbol l}=2.31'  with linespoints ls 3,\
'stereo_conv_l4.0000.dat'  using (fx($2)):(fy($4)) title '{/Symbol l}=4.00'  with linespoints ls 4,\
'stereo_conv_l5.0000.dat'  using (fx($2)):(fy($4)) title '{/Symbol l}=5.00'  with linespoints ls 5,\
'stereo_conv_l8.0000.dat'  using (fx($2)):(fy($4)) title '{/Symbol l}=8.00'  with linespoints ls 6,\
'stereo_conv_l16.0000.dat' using (fx($2)):(fy($4)) title '{/Symbol l}=16.00' with linespoints ls 7



set out    'G:\\LAT\\sd_metropolis\\plots\\rcan\\stereo_Gx_conv.eps'
set key right bottom
set xlabel "1/m_{max}"
set ylabel "<g_x>"
plot \
0 notitle with lines ls 0 lw 1,\
'stereo_conv_l0.5000.dat'  using (fx($2)):(fy($5)) title '{/Symbol l}=0.50'  with linespoints ls 2,\
'stereo_conv_l2.3100.dat'  using (fx($2)):(fy($5)) title '{/Symbol l}=2.31'  with linespoints ls 3,\
'stereo_conv_l4.0000.dat'  using (fx($2)):(fy($5)) title '{/Symbol l}=4.00'  with linespoints ls 4,\
'stereo_conv_l5.0000.dat'  using (fx($2)):(fy($5)) title '{/Symbol l}=5.00'  with linespoints ls 5,\
'stereo_conv_l8.0000.dat'  using (fx($2)):(fy($5)) title '{/Symbol l}=8.00'  with linespoints ls 6,\
'stereo_conv_l16.0000.dat' using (fx($2)):(fy($5)) title '{/Symbol l}=16.00' with linespoints ls 7
