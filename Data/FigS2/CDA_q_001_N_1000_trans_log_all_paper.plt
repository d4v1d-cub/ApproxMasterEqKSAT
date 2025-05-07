set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_q_001_N_1000_transition_log_all_paper.eps"
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


p for [i=1:8] 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2450_q_0.0010_tl_400.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 6 w l \
, for [i=1:8] 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2500_q_0.0010_tl_400.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 7 w l \
, for [i=1:8] 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2550_q_0.0010_tl_400.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 1 w l \
, for [i=1:8] 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2600_q_0.0010_tl_5000.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 2 w l \
, for [i=1:8] 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2700_q_0.0010_tl_5000.00_seed_'.i.'_tol_1.0e-03.txt' u 1:2 notitle lw 1.5 dt 2 lc 4 w l \
, 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2450_q_0.0010_tl_400.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.45}" lw 1.5 dt 2 lc 6 w l \
, 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2500_q_0.0010_tl_400.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.50}" lw 1.5 dt 2 lc 7 w l \
, 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2550_q_0.0010_tl_400.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.55}" lw 1.5 dt 2 lc 1 w l \
, 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2600_q_0.0010_tl_5000.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.60}" lw 1.5 dt 2 lc 2 w l \
, 'N_1000/CDA_WalkSAT_ener_K_3_N_1000_M_2700_q_0.0010_tl_5000.00_seed_1_tol_1.0e-03.txt' u 1:2 title "{/=24 {/Symbol a}=2.70}" lw 1.5 dt 2 lc 4 w l