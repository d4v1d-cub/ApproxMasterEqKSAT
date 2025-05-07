set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_vs_FMS_energy_eta_50_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set origin 0,0
set ylabel "{/=28 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=28 t}" offset 0,0
set label "{/=30 {/Symbol h} = 0.5}" at 250, 0.6
set xrange [0.01:2010]
set yrange [0.0001:1]
set key spacing 3.5 maxrows 10 width 6 at 1, 0.08
set tics  font ",24"
set ytics 1,.2,1
set ytics add ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "0.1" 0.1, "1" 1)


p 'FMS/FMS_KSAT_K_3_eta_50_N_50000_alpha_345_ngraphs_1_nhist_10000_tl_2000_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 5:5 title "{/=28 FMS {/Symbol a} = 3.45" ps 1.2 pt 4  lc 2 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_50_N_50000_alpha_358_ngraphs_1_nhist_10000_tl_2000_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 5:5 title "{/=28 {/Symbol a} = 3.58" ps 1.2 pt 6  lc 7 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_50_N_50000_alpha_365_ngraphs_1_nhist_10000_tl_2000_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 5:5 title "{/=28 {/Symbol a} = 3.65" ps 1.2 pt 8  lc 6 lw 2 w p \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_172500_eta_0.5000_tl_2000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 2 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_179000_eta_0.5000_tl_2000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 7 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_182500_eta_0.5000_tl_2000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 6 lw 3 w l \
, (200, 200) title "{/=28 CDA" lc 8 lw 2 w l
