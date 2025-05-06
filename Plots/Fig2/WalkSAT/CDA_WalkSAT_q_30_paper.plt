set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_vs_WalkSAT_energy_q_30_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set origin 0,0
set ylabel "{/=28 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=28 t}" offset 0,0
set label "{/=30 q = 0.3}" at 15, 0.6
set xrange [0.01:61]
set yrange [0.0001:1]
set key spacing 3.5 maxrows 10 width 6 at 0.2, 0.08
set tics  font ",24"
set ytics 1,.2,1
set ytics add ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "0.1" 0.1, "1" 1)


p 'WalkSAT/WalkSAT_dynamics_K_12_N_50000_alpha_2.55_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1<1?$1:0):2 every 10:10 title "{/=28 {/Symbol a} = 2.55" ps 1.5 pt 4  lc 2 lw 2 w p \
, 'WalkSAT/WalkSAT_dynamics_K_13_N_50000_alpha_2.7_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1<1?$1:0):2 every 10:10 title "{/=28 {/Symbol a} = 2.70" ps 1.6 pt 6  lc 7 lw 2 w p \
, 'WalkSAT/WalkSAT_dynamics_K_13_N_50000_alpha_2.92_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1<1?$1:0):2 every 10:10 title "{/=28 {/Symbol a} = 2.92" ps 1.5 pt 8  lc 6 lw 2 w p \
, 'WalkSAT/WalkSAT_dynamics_K_12_N_50000_alpha_2.55_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1>1?$1:0):2 every 1:1 notitle ps 1.5 pt 4  lc 2 lw 2 w p \
, 'WalkSAT/WalkSAT_dynamics_K_13_N_50000_alpha_2.7_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1>1?$1:0):2 every 1:1 notitle ps 1.6 pt 6  lc 7 lw 2 w p \
, 'WalkSAT/WalkSAT_dynamics_K_13_N_50000_alpha_2.92_q_0.3_t_60.0_dtmin_0.01_nhist_10000.txt'  u ($1>1?$1:0):2 every 2:2 notitle ps 1.5 pt 8  lc 6 lw 2 w p \
, 'CDA/CDA_WalkSAT_av_rates_ener_K_3_N_50000_M_127500_q_0.3000_tl_60.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 2 lw 3 w l \
, 'CDA/CDA_WalkSAT_av_rates_ener_K_3_N_50000_M_135000_q_0.3000_tl_60.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 7 lw 3 w l \
, 'CDA/CDA_WalkSAT_av_rates_ener_K_3_N_50000_M_146000_q_0.3000_tl_60.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 6 lw 3 w l \
, (200, 200) title "{/=28 CDA" lc 8 lw 2 w l
