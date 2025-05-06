set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_vs_FMS_energy_eta_70_paper.eps"
set output filenameoutput

reset
set size 1,1
set log
set origin 0,0
set ylabel "{/=28 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=28 t}" offset 0,0
set label "{/=30 {/Symbol h} = 0.7}" at 8, 0.6
set xrange [0.01:31]
set yrange [0.0001:1]
set key spacing 3.5 maxrows 10 width 6 at 0.2, 0.08
set tics  font ",24"
set ytics 1,.2,1
set ytics add ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "0.1" 0.1, "1" 1)


p 'FMS/FMS_KSAT_K_3_eta_70_N_50000_alpha_280_ngraphs_1_nhist_10000_tl_30_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 1:1 title "{/=28 {/Symbol a} = 2.80" ps 1.2 pt 4  lc 2 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_70_N_50000_alpha_300_ngraphs_1_nhist_10000_tl_30_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 1:1 title "{/=28 {/Symbol a} = 3.00" ps 1.2 pt 6  lc 7 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_70_N_50000_alpha_312_ngraphs_1_nhist_10000_tl_30_idgraph_1_logscale.txt'  u ($1/50000):($2/50000) every 1:1 title "{/=28 {/Symbol a} = 3.12" ps 1.2 pt 8  lc 6 lw 2 w p \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_140000_eta_0.7000_tl_30.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 2 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_150000_eta_0.7000_tl_30.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 7 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_50000_M_156000_eta_0.7000_tl_30.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 6 lw 3 w l \
, (200, 200) title "{/=28 CDA" lc 8 lw 2 w l
