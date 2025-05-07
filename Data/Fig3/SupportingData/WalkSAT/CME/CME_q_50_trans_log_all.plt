set term postscript enhanced color eps dl 2.5
filenameoutput="CME_q_50_transition_log_all.eps"
set output filenameoutput

reset
set multiplot
set size 1,1
set log
set ylabel "{/=24 e(t)}" rotate by 90 offset 1,3
set xlabel "{/=24 {t/N}}" 
set key spacing 2.0 maxrows 8 width 6 at 0.1,0.08
set label "{/=30 q = 0.5}" at 0.15, 0.6
set tics  font ",20"
set xrange[0.01:20.0]
set yrange[0.0001:1]
set ytics ("10^{-4}" 1e-4, "10^{-3}" 1e-3, "10^{-2}" 1e-2, "10^{-1}" 1e-1, 1)


do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_223000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 1 w l
}


do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_228000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 7 w l
}


do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_233000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 2 w l
}


do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_234000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 3 w l
}

do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_235000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 4 w l
}


do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_236000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 8 w l
}

do for [i=1:8]{
filename = sprintf('N_100000/CME_WalkSAT_av_rates_ener_K_3_N_100000_M_238000_q_0.5000_tl_20.00_seed_%d_tol_1.0e-03.txt',i)
p filename u 1:2 notitle lw 1.5 dt 2 lc 6 w l
}