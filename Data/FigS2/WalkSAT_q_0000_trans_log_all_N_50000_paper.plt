set term postscript enhanced color eps dl 2.5
filenameoutput="WalkSAT_q_0000_transition_log_all_N_50000_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.07,0.0002
set label "{/=30 q = 0.0}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:10.0]
set yrange[0.000001:1]
set ytics ("10^{-6}" 1e-6, "10^{-5}" 1e-5, "10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, "1"1)

p 'N_50000/WalkSAT_dynamics_N_50000_alpha_0.5_q_0.0_t_5000.0_dtmin_0.01.txt' u 1:2 title "{/=24 {/Symbol a}=0.5}" lw 2 dt 1 lc 4 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_0.6_q_0.0_t_20.0_dtmin_0.01_ngraphs_1000_nhist_1.txt' u 1:2 title "{/=24 {/Symbol a}=0.6}" lw 2 dt 1 lc 6 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_0.8_q_0.0_t_20.0_dtmin_0.01_ngraphs_1000_nhist_1.txt' u 1:2 title "{/=24 {/Symbol a}=0.8}" lw 2 dt 1 lc 7 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_1.0_q_0.0_t_20.0_dtmin_0.01_ngraphs_1000_nhist_1.txt' u 1:2 title "{/=24 {/Symbol a}=1.0}" lw 2 dt 1 lc 8 w l \
, 'N_50000/WalkSAT_dynamics_N_50000_alpha_1.2_q_0.0_t_20.0_dtmin_0.01_ngraphs_1000_nhist_1.txt' u 1:2 title "{/=24 {/Symbol a}=1.2}" lw 2 dt 1 lc 1 w l