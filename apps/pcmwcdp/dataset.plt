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

eval('log_link_f1(x) = log_link_'.msuffix.'_extfunc(x) - sqrt(log_link_'.msuffix.'_cfdfunc(x))')
eval('log_link_f2(x) = log_link_'.msuffix.'_extfunc(x) + sqrt(log_link_'.msuffix.'_cfdfunc(x))')
eval('log_link_ext   = log_link_'.msuffix.'_ext')
eval('log_link_err   = log_link_'.msuffix.'_err')

eval('log_mpnt_f1(x) = log_mpnt_'.msuffix.'_extfunc(x) - sqrt(log_mpnt_'.msuffix.'_cfdfunc(x))')
eval('log_mpnt_f2(x) = log_mpnt_'.msuffix.'_extfunc(x) + sqrt(log_mpnt_'.msuffix.'_cfdfunc(x))')
eval('log_mpnt_ext   = log_mpnt_'.msuffix.'_ext')
eval('log_mpnt_err   = log_mpnt_'.msuffix.'_err')

eval('log_Gx_f1(x)   = log_Gx_'.msuffix.'_extfunc(x) - sqrt(log_Gx_'.msuffix.'_cfdfunc(x))')
eval('log_Gx_f2(x)   = log_Gx_'.msuffix.'_extfunc(x) + sqrt(log_Gx_'.msuffix.'_cfdfunc(x))')
eval('log_Gx_ext     = log_Gx_'.msuffix.'_ext')
eval('log_Gx_err     = log_Gx_'.msuffix.'_err')

eval('srt_link_f1(x) = sqrt_link_'.msuffix.'_extfunc(x) - sqrt(sqrt_link_'.msuffix.'_cfdfunc(x))')
eval('srt_link_f2(x) = sqrt_link_'.msuffix.'_extfunc(x) + sqrt(sqrt_link_'.msuffix.'_cfdfunc(x))')
eval('srt_link_ext   = sqrt_link_'.msuffix.'_ext')
eval('srt_link_err   = sqrt_link_'.msuffix.'_err')

eval('srt_mpnt_f1(x) = sqrt_mpnt_'.msuffix.'_extfunc(x) - sqrt(sqrt_mpnt_'.msuffix.'_cfdfunc(x))')
eval('srt_mpnt_f2(x) = sqrt_mpnt_'.msuffix.'_extfunc(x) + sqrt(sqrt_mpnt_'.msuffix.'_cfdfunc(x))')
eval('srt_mpnt_ext   = sqrt_mpnt_'.msuffix.'_ext')
eval('srt_mpnt_err   = sqrt_mpnt_'.msuffix.'_err')

eval('srt_Gx_f1(x)   = sqrt_Gx_'.msuffix.'_extfunc(x) - sqrt(sqrt_Gx_'.msuffix.'_cfdfunc(x))')
eval('srt_Gx_f2(x)   = sqrt_Gx_'.msuffix.'_extfunc(x) + sqrt(sqrt_Gx_'.msuffix.'_cfdfunc(x))')
eval('srt_Gx_ext     = sqrt_Gx_'.msuffix.'_ext')
eval('srt_Gx_err     = sqrt_Gx_'.msuffix.'_err')

plot_dir = 'G:\\LAT\\sd_metropolis\\plots\\pcm_wc_mspace'

set xlabel "1/M"
set yrange [-0.01:]
set xrange [0.0:1.02]
set key right bottom

set out sprintf("%s\\Gx_convergence_%s.eps", plot_dir, suffix)
set ylabel "<1/N tr g_x>_M"
plot \
'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact' using ($1):($2)      notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'Gx_'.suffix.'.mean'                                using ($1):($2):($3) notitle with yerrorbars ls 2,\
'+' using     1:(log_Gx_f1($1)):(log_Gx_f2($1)) notitle          with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(log_Gx_ext)   :(log_Gx_err)    title 'Log ext'  with yerrorbars   ls 2 ps 3,\
'+' using     1:(srt_Gx_f1($1)):(srt_Gx_f2($1)) notitle          with filledcurves ls 2 lc rgb '#AAAAFF',\
'+' using (0.0):(srt_Gx_ext)   :(srt_Gx_err)    title 'Sqrt ext' with yerrorbars   ls 4 ps 3

set out sprintf("%s\\mpnt_convergence_%s.eps", plot_dir, suffix)
set ylabel "G_{xy}(L_s/2)"
plot \
'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact'   using ($1):($4)                     notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc6'  using (0.0):($4) every ::54::54     title 'N=6' with points     ls 2 lc rgb '#AAFFAA' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc9'  using (0.0):($4) every ::54::54     title 'N=9' with points     ls 3 lc rgb '#77FF77' ps 3,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc15' using (0.0):($4) every ::54::54     title 'N=15' with points    ls 4 lc rgb '#22FF22' ps 3,\
'mpnt_'.suffix.'.mean'                                using ($1):($2):($3)                notitle with yerrorbars ls 2,\
'+' using     1:(log_mpnt_f1($1)):(log_mpnt_f2($1))  notitle            with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(log_mpnt_ext)   :(log_mpnt_err)     title 'Log ext'    with yerrorbars ls 2 ps 3,\
'+' using     1:(srt_mpnt_f1($1)):(srt_mpnt_f2($1))  notitle            with filledcurves ls 2 lc rgb '#AAAAFF',\
'+' using (0.0):(srt_mpnt_ext)   :(srt_mpnt_err)     title 'Sqrt ext'   with yerrorbars ls 4 ps 3

#Three-loop result for mean link, taken from hep-lat/9401029
a1 = 0.03125
a2 = 0.00435762
E1 = 1.0 - 0.125*lambda
E2 = 1.0 - 0.125*lambda*(1.0 + a1*lambda)
E3 = 1.0 - 0.125*lambda*(1.0 + a1*lambda + a2*lambda*lambda)

xnc(n) = 10.0/(n*n)

#Infinite-N extrapolations for publication plots
if(abs(lambda-3.0120)<0.001) { sun_link_ext(x) = (0.544253 + 0.80996*x) } else {}
if(abs(lambda-3.1000)<0.001) { sun_link_ext(x) = (0.516138 + 1.18254*x) } else {}
if(abs(lambda-3.2300)<0.001) { sun_link_ext(x) = (0.482232 + 1.31012*x) } else {}
if(abs(lambda-4.0000)<0.001) { sun_link_ext(x) = (0.284807 + 1.15307*x) } else {}

set arrow from 0.0,sun_link_ext(0.0) to 10.0/(6*6),sun_link_ext(1.0/(6*6)) nohead ls 3 lw 2

set autoscale y
set out sprintf("%s\\link_convergence_%s.eps", plot_dir, suffix)
set ylabel "<1/N tr(g_1^+g_0)>"
set xlabel "1/M,   10/N^2"
set label sprintf("{/Symbol l} = %2.4f", lambda) at graph 0.5,0.9 center
plot \
1/0 title 'MC,SU(N)' with linespoints ls 3 ps 2.0 lw 2,\
for [nc = 6:21:3] \
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc'.nc  using  (xnc(nc)):($4) every ::1::1 notitle  with points     ls 3 ps 2.0,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.unc9'   using  (xnc(9)):($4) every ::1::1  title 'MC,U(N)'  with points     ls 5 ps 2.0,\
for [nc = 15:21:6] \
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.unc'.nc  using  (xnc(nc)):($4) every ::1::1 notitle  with points     ls 5 ps 2.0,\
'+' using (1.00):(E1) title 'Standard PT' with points ls 1 ps 1.7,\
'+' using (0.50):(E2) notitle             with points ls 1 ps 1.7,\
'+' using (0.33):(E3) notitle             with points ls 1 ps 1.7,\
'link_'.suffix.'.mean'                              using ($1):($2):($3) title 'DiagMC' with yerrorbars ls 2,\
'+' using     1:(log_link_f1($1)):(log_link_f2($1)) notitle          with filledcurves ls 2 lc rgb '#FFAAAA',\
'+' using (0.0):(log_link_ext)   :(log_link_err)    title 'Log ext'  with yerrorbars   ls 2 ps 1.8,\
'+' using     1:(srt_link_f1($1)):(srt_link_f2($1)) notitle          with filledcurves ls 2 lc rgb '#AAAAFF',\
'+' using (0.0):(srt_link_ext)   :(srt_link_err)    title 'Sqrt ext' with yerrorbars   ls 4 ps 1.8

unset  label
unset  arrow

#'..\\pcm_wc_mspace_exact\\scalars_'.suffix.'.exact' using ($1):($3)      notitle with points     ls 2 lc rgb '#AAFFAA' ps 3,\

set xlabel "1/M"

#set out sprintf("%s\\correlator_%s.eps", plot_dir, suffix)
#set yrange [-0.01:1.01]
#set xrange [-0.2:20.2]
#set xlabel "x"
#set ylabel "<1/N tr(g_x^+g_0)>_M"
#plot \
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.exact' using ($1):($2) notitle with lines lt 1 lc rgbcolor word(blueredcolortable, 1)  lw 4,\
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.exact' using ($1):($3) notitle with lines lt 1 lc rgbcolor word(blueredcolortable, 2)  lw 4,\
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc6'   using ($3):($4) title 'N=6'  with lines lt 1 lc rgb '#00AA00' lw 4,\
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc9'   using ($3):($4) title 'N=9'  with lines lt 1 lc rgb '#00DD00' lw 4,\
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc12'  using ($3):($4) title 'N=12' with lines lt 1 lc rgb '#00EE00' lw 4,\
#'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc18'  using ($3):($4) title 'N=18' with lines lt 1 lc rgb '#00FF00' lw 4,\
#for [io=0:11] \
#'Gxy_'.suffix.'.dat' using ($1):(column(2*io+2)):(column(2*io+3)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7

#Analytic results for Green's functions, hep-lat/9401029
EVR(x) = (0.5/pi)*(-log(x) - 0.57721 - 1.5*log(2.0))
G1(x) = 1 + 0.5*lambda*EVR(x)
G2(x) = G1(x) + (lambda*lambda/32.0)*EVR(x)*(1.0 + 2.0*EVR(x))

nc1 = 0
nc2 = 0
if(abs(lambda-3.0120)<0.001) { nc1 = 12; nc2 = 18; } else {}
if(abs(lambda-3.1000)<0.001) { nc1 = 6; nc2 = 9; } else {}

set out sprintf("%s\\correlator_%s.eps", plot_dir, suffix)
set label sprintf("{/Symbol l} = %2.4f", lambda) at graph 0.2,0.1 center
set key top right maxrows 3
set yrange [-0.2:1.1]
set xrange [-0.2:25.2]
set xlabel "x"
set ylabel "<1/N tr(g_x^+g_0)>"
plot \
'+' using (0.0):(2.0):(0.0) title 'DiagMC,M=1'  with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, 1)  lw 4 pt 7,\
'+' using (0.0):(2.0):(0.0) title 'DiagMC,M=12' with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, 12) lw 4 pt 7,\
0 notitle with lines ls 7 lw 0.8,\
G1(x) title 'Std.PT,LO'  with lines ls 1,\
G2(x) title 'Std.PT,NLO' with lines ls 7,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc'.nc1  using ($3):($4):($5) title 'MC,SU('.nc1.')' with yerrorbars ls 3 ps 1.7,\
'..\\pcm_wc_mspace_exact\\correlator_'.suffix.'.nc'.nc2  using ($3):($4):($5) title 'MC,SU('.nc2.')' with yerrorbars ls 5 ps 1.7,\
for [io=0:11] \
'Gxy_'.suffix.'.dat' using ($1):(column(2*io+2)):(column(2*io+3)) notitle with yerrorbars lt 1 lc rgbcolor word(blueredcolortable, io+1) lw 4 pt 7

unset label

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


