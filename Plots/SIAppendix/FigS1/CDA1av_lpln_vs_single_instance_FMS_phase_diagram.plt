set term postscript enhanced color eps dl 2.5
filenameoutput="CDA1av_lpln_vs_single_instance_FMS_phase_diagram.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=26 {/Symbol h}}" rotate by 90 offset -1,-1
set xlabel "{/=26 {/Symbol a}}" offset 0,0
set key spacing 2.0 maxrows 4 width 6 at 2.27,0.4
set xrange [1.0:5.5]
set yrange [0.02:1.01]
set tics  font ",24"
set arrow nohead from 3.86,0.02 to 3.86,1.01 lc 8 lw 2
set arrow nohead from 4.27,0.02 to 4.27,1.01 lc 8 lw 2
set label "{/=28 SAT" at 3.5, 0.15
set label "{/=28 UNSAT" at 4.35, 0.15
set label "{/=26 {/Symbol a}_d" at 3.75, -0.03
set label "{/=26 {/Symbol a}_s" at 4.2, -0.03
set ytics 0, 0.1, 1

p 'FMS_KSAT_K_3_PD_convexity.txt' u 2:1:($2-$3):($2+$4) w xerrorbars ps 1.9 pt 8 lc 8 lw 2 title "{/=26 FMS" \
, 'CDA_transition_FMS_N_1000.txt' u 2:1:($2-$3):($2+$4) w xerrorbars ps 1 pt 4 lc 2 lw 3 title "{/=26 CDA" \
, 'CDA_transition_FMS_N_1000.txt' u 2:1 w l dt 2 lc 2 lw 2 notitle \
, 'CME_KSAT_K_3_phase_diagram_N_10000_newcode.txt' u 2:1:($2-$3):($2+$4) w xerrorbars ps 1.1 pt 6 lc 6 lw 3 title "{/=28 CME" \
, 'CME_KSAT_K_3_phase_diagram_N_10000_newcode.txt' u 2:1 w l dt 2 lc 6 lw 2 notitle \
, 'CDA1av_lpln_popdyn_transitions.txt' u 2:1:($2-$3):($2+$4) w xerrorbars ps 1 pt 1 lc 2 lw 4 title "{/=26 av-CDA" \
, 'CDA1av_lpln_popdyn_transitions.txt' u 2:1 w l lc 2 lw 2 dt 1 notitle
