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

load 'G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace\\extrapolations.plt'

eval('link_f1(x) = link_'.msuffix.'_extfunc(x) - sqrt(link_'.msuffix.'_cfdfunc(x))')
eval('link_f2(x) = link_'.msuffix.'_extfunc(x) + sqrt(link_'.msuffix.'_cfdfunc(x))')
eval('link_ext = link_'.msuffix.'_ext')
eval('link_err = link_'.msuffix.'_err')

eval('mpnt_f1(x) = mpnt_'.msuffix.'_extfunc(x) - sqrt(mpnt_'.msuffix.'_cfdfunc(x))')
eval('mpnt_f2(x) = mpnt_'.msuffix.'_extfunc(x) + sqrt(mpnt_'.msuffix.'_cfdfunc(x))')
eval('mpnt_ext = mpnt_'.msuffix.'_ext')
eval('mpnt_err = mpnt_'.msuffix.'_err')

eval('Gx_f1(x) = Gx_'.msuffix.'_extfunc(x) - sqrt(Gx_'.msuffix.'_cfdfunc(x))')
eval('Gx_f2(x) = Gx_'.msuffix.'_extfunc(x) + sqrt(Gx_'.msuffix.'_cfdfunc(x))')
eval('Gx_ext = Gx_'.msuffix.'_ext')
eval('Gx_err = Gx_'.msuffix.'_err')

plot_dir = 'G:\\LAT\\sd_metropolis\\plots\\pcm_wc_mspace'

set xlabel "1/m_{max}"
set yrange [-0.01:]
set xrange [0.0:1.02]
set key right bottom

set out sprintf("%s\\Gx_convergence_%s.eps", plot_dir, suffix)
set ylabel "G_x"
plot \
'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact' using ($1):($2)      notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'Gx_'.suffix.'.mean'                                using ($1):($2):($3) notitle with yerrorbars ls 2,\
'+' using 1:(Gx_f1($1)):(Gx_f2($1)) notitle with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(Gx_ext):(Gx_err) notitle with yerrorbars ls 4 ps 3

set out sprintf("%s\\mpnt_convergence_%s.eps", plot_dir, suffix)
set ylabel "G_{xy}(L_s/2)"
plot \
'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact'   using ($1):($4)                     notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc6'  using (0.0):($4) every ::54::54     title 'N=6' with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc9'  using (0.0):($4) every ::54::54     title 'N=9' with points     ls 3 lc rgb '#77FF77' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc15' using (0.0):($4) every ::54::54     title 'N=15' with points    ls 4 lc rgb '#22FF22' ps 3,\
'mpnt_'.suffix.'.mean'                                using ($1):($2):($3)                notitle with yerrorbars ls 2,\
'+' using 1:(mpnt_f1($1)):(mpnt_f2($1)) notitle with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(mpnt_ext):(mpnt_err) title 'Ext' with yerrorbars ls 4 ps 3


set autoscale y
set out sprintf("%s\\link_convergence_%s.eps", plot_dir, suffix)
set ylabel "Link"
plot \
'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact' using ($1):($3)      notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc6'  using (0.0):($4) every ::1::1     title 'N=6'  with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc9'  using (0.0):($4) every ::1::1     title 'N=9'  with points     ls 3 lc rgb '#77FF77' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc15' using (0.0):($4) every ::1::1     title 'N=15' with points     ls 4 lc rgb '#22FF22' ps 3,\
'link_'.suffix.'.mean'                              using ($1):($2):($3) notitle with yerrorbars ls 2,\
'+' using 1:(link_f1($1)):(link_f2($1)) notitle with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(link_ext):(link_err) title 'Ext' with yerrorbars ls 4 ps 3


set out sprintf("%s\\correlator_%s.eps", plot_dir, suffix)
set yrange [-0.01:1.01]
set xrange [-0.2:54.2]
set xlabel "x"
set ylabel "G_{xy}(x)"
plot \
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.exact' using ($1):($2) notitle with lines lt 1 lc rgbcolor word(blueredcolortable, 1)  lw 4,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.exact' using ($1):($3) notitle with lines lt 1 lc rgbcolor word(blueredcolortable, 2)  lw 4,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc6'   using ($3):($4) title 'N=6' with lines lt 1 lc rgb '#00AA00' lw 4,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc9'   using ($3):($4) title 'N=9' with lines lt 1 lc rgb '#00DD00' lw 4,\
for [io=0:11] \
'Gxy_'.suffix.'.dat' using ($1):(column(2*io+2)):(column(2*io+3)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7

set autoscale y
set logscale y
set autoscale x
set key left bottom
set out sprintf("%s\\signs_%s.eps", plot_dir, suffix)
set ylabel "Mean sign"
set xlabel "Expansion order"
plot \
'signs_'.suffix.'.mean' using ($1):(abs($2)):($3) title 'Gx' with yerrorbars   ls 2,\
'signs_'.suffix.'.mean' using ($1):(abs($4)):($5) title 'Link' with yerrorbars ls 3

quit


