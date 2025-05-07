set term postscript enhanced color eps dl 2.5
filenameoutput="WalkSAT_q_001_transition_log_all_N_10000_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.3,0.008
set label "{/=30 q = 0.001}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:500.0]
set yrange[0.0001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, "1"1)

p 'N_10000/WalkSAT_dynamics_N_10000_alpha_1.5_q_0.001_t_1000.0.txt' u 1:2 ev 50:50 title "{/=24 {/Symbol a}=1.5}" lw 2 dt 1 lc 1 w l \
, 'N_10000/WalkSAT_dynamics_N_10000_alpha_1.9_q_0.001_t_1000.0.txt' u 1:2 ev 50:50 title "{/=24 {/Symbol a}=1.9}" lw 2 dt 1 lc 2 w l \
, 'N_10000/WalkSAT_dynamics_N_10000_alpha_2.3_q_0.001_t_1000.0.txt' u 1:2 ev 50:10 title "{/=24 {/Symbol a}=2.3}" lw 2 dt 1 lc 4 w l \
, 'N_10000/WalkSAT_dynamics_N_10000_alpha_2.7_q_0.001_t_1000.0.txt' u 1:2 ev 50:50 title "{/=24 {/Symbol a}=2.7}" lw 2 dt 1 lc 6 w l \
, 'N_10000/WalkSAT_dynamics_N_10000_alpha_2.9_q_0.001_t_1000.0.txt' u 1:2 ev 50:50 title "{/=24 {/Symbol a}=2.9}" lw 2 dt 1 lc 7 w l