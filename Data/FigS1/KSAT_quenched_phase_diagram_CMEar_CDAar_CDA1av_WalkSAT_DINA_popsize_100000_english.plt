set term postscript enhanced color eps dl 2.5
filenameoutput="CMEar_CDAar_CDA1av_DINA_WalkSAT_KSAT_popsize_100000_quenched_phase_diagram_english.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=26 q}" rotate by 90 offset -1,-1
set xlabel "{/=26 {/Symbol a}}" offset 0,0
set key spacing 2.0 maxrows 5 width 6 at 2.27,0.4
set xrange [1.0:4.6]
set yrange [0.02:1.02]
set tics  font ",24"
set arrow nohead from 3.86,0.02 to 3.86,1.02 lc 8 lw 2
set arrow nohead from 4.27,0.02 to 4.27,1.02 lc 8 lw 2
set label "{/=26 G-WalkSAT" at 1.21, 0.85
set label "{/=26 succeeds" at 1.275, 0.79
set label "{/=26 G-WalkSAT" at 2.962, 0.85
set label "{/=26 fails" at 3.22, 0.79
set label "{/=26 {/Symbol a}_d" at 3.75, -0.03
set label "{/=26 {/Symbol a}_s" at 4.2, -0.03
set ytics 0, 0.1, 1

p 'WalkSAT_PD_convexity.txt' u (($2 + $3)/2):1:2:3 w xerror ps 1.9 pt 8 lc 8 lw 3  title "{/=26 G-WalkSAT" \
, 'RE_N-10000 _K-3 _size-25 diagrama de fases_mano' u (($2 + $3)/2):1:2:3 w xerror lc 7 lw 2.5 title "{/=26 DINA" \
, 'CDA_av_rates_K_3_PD.txt' u (($2 + $3)/2):1:2:3 w xerror ps 1 pt 4 lc 2 lw 3 title "{/=26 CDA" \
, 'CME_av_rates_K_3_PD.txt' u (($2 + $3)/2):1:2:3 w xerror ps 1.1 pt 6 lc 6 lw 3 title "{/=26 CME" \
, 'CDA1av_lpln_popsize_100000_K_3_PD.txt' u (($2 + $3)/2):1:2:3 w xerror ps 1 pt 1 lc 2 lw 3 title "{/=26 av-CDA" \
, 'RE_N-10000 _K-3 _size-25 diagrama de fases_mano' u (($2 + $3)/2):1 w l lc 7 lw 1.5 dt 2 notitle \
, 'CDA_av_rates_K_3_PD.txt' u (($2 + $3)/2):1 w l lc 2 lw 2 dt 2 notitle \
, 'CME_av_rates_K_3_PD.txt' u (($2 + $3)/2):1 w l lc 6 lw 2 dt 2 notitle \
, 'CDA1av_lpln_popsize_100000_K_3_PD.txt' u (($2 + $3)/2):1 w l dt 1 lc 2 lw 2 notitle