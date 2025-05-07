set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_vs_FMS_energy_eta_30_N_500000.eps"
set output filenameoutput

reset
set size 1,1
set log
set origin 0,0
set ylabel "{/=28 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=28 t}" offset 0,0
set label "{/=30 {/Symbol h} = 0.3}" at 250, 0.6
set xrange [0.01:2010]
set yrange [0.0001:1]
set key spacing 3.5 maxrows 10 width 6 at 1, 0.08
set tics  font ",24"
set ytics 1,.2,1
set ytics add ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "0.1" 0.1, "1" 1)


p 'FMS/FMS_KSAT_K_3_eta_30_N_500000_alpha_405_ngraphs_1_nhist_10000_tl_10000_idgraph_1_logscale.txt'  u ($1/500000):($2/500000) every 1:1 title "{/=28 FMS {/Symbol a} = 4.05" ps 1.2 pt 3  lc 1 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_30_N_500000_alpha_410_ngraphs_1_nhist_10000_tl_10000_idgraph_1_logscale.txt'  u ($1/500000):($2/500000) every 1:1 title "{/=28 {/Symbol a} = 4.10" ps 1.2 pt 4  lc 2 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_30_N_100000_alpha_415_ngraphs_1_nhist_10000_tl_10000_idgraph_1_logscale.txt'  u ($1/100000):($2/100000) every 1:1 title "{/=28 {/Symbol a} = 4.15" ps 1.2 pt 8  lc 6 lw 2 w p \
, 'FMS/FMS_KSAT_K_3_eta_30_N_500000_alpha_420_ngraphs_1_nhist_10000_tl_10000_idgraph_1_logscale.txt'  u ($1/500000):($2/500000) every 1:1 title "{/=28 {/Symbol a} = 4.20" ps 1.2 pt 10  lc 7 lw 2 w p \
, 'CDA/CDA_FMS_ener_K_3_N_500000_M_2025000_eta_0.3000_tl_10000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 1 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_500000_M_2050000_eta_0.3000_tl_10000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 2 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_500000_M_2075000_eta_0.3000_tl_10000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 6 lw 3 w l \
, 'CDA/CDA_FMS_ener_K_3_N_500000_M_2100000_eta_0.3000_tl_10000.00_seed_1_tol_1.0e-03.txt'  u 1:2  notitle lc 7 lw 3 w l \
, (200, 200) title "{/=28 CDA" lc 8 lw 2 w l
