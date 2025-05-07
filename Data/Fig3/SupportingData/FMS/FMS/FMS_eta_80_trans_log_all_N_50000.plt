set term postscript enhanced color eps dl 2.5
filenameoutput="FMS_eta_80_transition_log_all_N_50000.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 {/Symbol h} = 0.8}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.1:110.0]
set yrange[0.0001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, 1)

p 'N_50000/FMS_dynamics_K_3_N_50000_eta_80_alpha_290_tl_100_dtN_10000_nsamples_5000.txt' u 1:($2/50000) notitle lw 1.5 dt 1 lc 1 w l \
, 'N_50000/FMS_dynamics_K_3_N_50000_eta_80_alpha_292_tl_100_dtN_10000_nsamples_3535.txt' u 1:($2/50000) notitle lw 1.5 dt 1 lc 2 w l \
, 'N_50000/FMS_dynamics_K_3_N_50000_eta_80_alpha_294_tl_50_dtN_5000_nsamples_1000.txt' u 1:($2/50000) notitle lw 1.5 dt 1 lc 6 w l \
, 'N_50000/FMS_dynamics_K_3_N_50000_eta_80_alpha_296_tl_50_dtN_5000_nsamples_1000.txt' u 1:($2/50000) notitle lw 1.5 dt 1 lc 7 w l