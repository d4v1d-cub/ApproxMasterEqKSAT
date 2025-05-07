set term postscript enhanced color eps dl 2.5
filenameoutput="FMS_eta_10_transition_log_all_N_100000.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 {/Symbol h} = 0.1}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.1:550000.0]
set yrange[0.000001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, 1)

p 'N_100000/FMS_dynamics_K_3_N_100000_eta_10_alpha_400_tl_100000_dtN_10000000_nsamples_1000.txt' u 1:($2/100000) notitle lw 1.5 dt 1 lc 1 w l \
, 'N_100000/FMS_dynamics_K_3_N_100000_eta_10_alpha_402_tl_500000_dtN_10000000_nsamples_1000.txt' u 1:($2/100000) notitle lw 1.5 dt 1 lc 2 w l \
, 'N_100000/FMS_dynamics_K_3_N_100000_eta_10_alpha_404_tl_500000_dtN_10000000_nsamples_166.txt' u 1:($2/100000) notitle lw 1.5 dt 1 lc 6 w l \
, 'N_100000/FMS_dynamics_K_3_N_100000_eta_10_alpha_406_tl_500000_dtN_10000000_nsamples_61.txt' u 1:($2/100000) notitle lw 1.5 dt 1 lc 7 w l