set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_q_0000_N_5000_transition_log_all_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.08,0.008
set label "{/=30 q = 0.0}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:10.0]
set yrange[0.000001:1]
set ytics ("10^{-6}" 1e-6, "10^{-5}" 1e-5, "10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, "1"1)


p for [i=1:8] 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_6000_q_0.0000_tl_10.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 6 w l \
, for [i=1:8] 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_8000_q_0.0000_tl_10.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 7 w l \
, for [i=1:8] 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_10000_q_0.0000_tl_10.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 1 w l \
, for [i=1:8] 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_12000_q_0.0000_tl_10.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 2 w l \
, for [i=1:8] 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_13000_q_0.0000_tl_10.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 4 w l \
, 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_6000_q_0.0000_tl_10.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=1.2}" lw 1.5 dt 2 lc 6 w l \
, 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_8000_q_0.0000_tl_10.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=1.6}" lw 1.5 dt 2 lc 7 w l \
, 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_10000_q_0.0000_tl_10.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.0}" lw 1.5 dt 2 lc 1 w l \
, 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_12000_q_0.0000_tl_10.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.4}" lw 1.5 dt 2 lc 2 w l \
, 'N_5000/CDA_WalkSAT_ener_K_3_N_5000_M_13000_q_0.0000_tl_10.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.6}" lw 1.5 dt 2 lc 4 w l