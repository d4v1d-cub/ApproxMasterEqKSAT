set term postscript enhanced color eps dl 2.5
filenameoutput="WalkSAT_q_40_transition_log_all_N_10000.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 q = 0.4}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:20.0]
set yrange[0.0001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, 1)

p 'N_10000/WalkSAT_dynamics_N_10000_alpha_2.87_q_0.4_t_10.0.txt' u 1:2 notitle lw 1.5 dt 1 lc 1 w l \
, 'N_10000/WalkSAT_dynamics_N_10000_alpha_2.9_q_0.4_t_10.0.txt' u 1:2 notitle lw 1.5 dt 1 lc 2 w l