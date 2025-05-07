set term postscript enhanced color eps dl 2.5
filenameoutput="WalkSAT_q_90_transition_log_all_N_50000.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 q = 0.9}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:60.0]
set yrange[0.0001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, 1)

p 'N_50000/WalkSAT_dynamics_N_50000_alpha_2.7_q_0.9_t_50.0_dtmin_0.01.txt' u 1:2 notitle lw 1.5 dt 1 lc 1 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_2.73_q_0.9_t_50.0_dtmin_0.01.txt' u 1:2 notitle lw 1.5 dt 1 lc 6 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_2.75_q_0.9_t_50.0_dtmin_0.01.txt' u 1:2 notitle lw 1.5 dt 1 lc 2 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_2.8_q_0.9_t_50.0_dtmin_0.01.txt' u 1:2 notitle lw 1.5 dt 1 lc 4 w l