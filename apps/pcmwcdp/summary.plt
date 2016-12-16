reset
set term postscript color enhanced landscape "Helvetica" 24
cd 'G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace\\'
set pointsize 1.0
set bar 1.0
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

blueredcolortable = "#0000ff #1700e8 #2e00d1 #4600b9 #5d00a2 #74008b #8b0074 #a2005d #b90046 #d1002e #e80017 #ff0000"

plot_dir = 'G:\\LAT\\sd_metropolis\\plots\\pcm_wc_mspace'

set autoscale x
set autoscale y
set xlabel suffix

set key right bottom
set out sprintf("%s\\acceptance_summary_%s.eps", plot_dir, suffix)
set ylabel "Acceptance"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($2):($3) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\mean_recursion_depth_summary_%s.eps", plot_dir, suffix)
set ylabel "Mean rec.depth"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($4):($5) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\mean_stack_top_summary_%s.eps", plot_dir, suffix)
set ylabel "Mean stack top"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($8):($9) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\mean_nA_summary_%s.eps", plot_dir, suffix)
set ylabel "1 - <N_A>"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($12):($13) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\mean_sign_summary_%s.eps", plot_dir, suffix)
set ylabel "Mean sign"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($16):($17) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\mean_return_time_summary_%s.eps", plot_dir, suffix)
set ylabel "Mean return time"
plot \
'mc_stat_summary_'.suffix.'.dat' using ($1):($18):($19) notitle with yerrorbars ls 2

set key right bottom
set out sprintf("%s\\Gx_summary_%s.eps", plot_dir, suffix)
set ylabel "Gx"
plot \
for [io=0:11] \
'Gx_summary_'.suffix.'.dat' using ($1):(column(2+2*io)):(column(3+2*io)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7,\
'Gx_summary_'.suffix.'.ext' using ($1):($2):($3) notitle with yerrorbars ls 5

set key right bottom
set out sprintf("%s\\link_summary_%s.eps", plot_dir, suffix)
set ylabel "Mean link"
plot \
'..\\pcm_wc_mspace_exact\\link_summary_'.suffix.'.nc6'  using ($1):($2) notitle with points ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\link_summary_'.suffix.'.nc9'  using ($1):($2) notitle with points ls 2 lc rgb '#77FF77' ps 3,\
'..\\pcm_wc_mspace_exact\\link_summary_'.suffix.'.nc15' using ($1):($2) notitle with points ls 2 lc rgb '#22FF22' ps 3,\
for [io=0:11] \
'link_summary_'.suffix.'.dat'                           using ($1):(column(2+2*io)):(column(3+2*io)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7,\
'link_summary_'.suffix.'.ext' using ($1):($2):($3) notitle with yerrorbars ls 5

set key right bottom
set out sprintf("%s\\mpnt_summary_%s.eps", plot_dir, suffix)
set ylabel "Correlator midpoint"
plot \
'..\\pcm_wc_mspace_exact\\mpnt_summary_'.suffix.'.nc6'  using ($1):($2) notitle with points ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\mpnt_summary_'.suffix.'.nc9'  using ($1):($2) notitle with points ls 2 lc rgb '#77FF77' ps 3,\
'..\\pcm_wc_mspace_exact\\mpnt_summary_'.suffix.'.nc15' using ($1):($2) notitle with points ls 2 lc rgb '#22FF22' ps 3,\
for [io=0:11] \
'mpnt_summary_'.suffix.'.dat'                           using ($1):(column(2+2*io)):(column(3+2*io)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7,\
'mpnt_summary_'.suffix.'.ext' using ($1):($2):($3) notitle with yerrorbars ls 5

quit
