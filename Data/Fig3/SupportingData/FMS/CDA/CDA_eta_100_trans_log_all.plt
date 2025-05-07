set term postscript enhanced color eps dl 2.5
filenameoutput="CDA_eta_100_transition_log_all.eps"
set output filenameoutput

reset
set multiplot
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 {/Symbol h} = 1.0}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:30.0]
set yrange[0.000001:1]
set ytics ("10^{-6}" 1e-6, "10^{-5}" 1e-5, "10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, "1" 1)


do for [i=1:8]{
filename = sprintf('N_10000/CDA_FMS_ener_K_3_N_10000_M_24600_eta_1.0000_tl_30.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 1 w l
}


do for [i=1:8]{
filename = sprintf('N_10000/CDA_FMS_ener_K_3_N_10000_M_24900_eta_1.0000_tl_30.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 2 w l
}


do for [i=1:8]{
filename = sprintf('N_10000/CDA_FMS_ener_K_3_N_10000_M_25100_eta_1.0000_tl_30.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 6 w l
}


do for [i=1:8]{
filename = sprintf('N_10000/CDA_FMS_ener_K_3_N_10000_M_25400_eta_1.0000_tl_30.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 7 w l
}