set term postscript enhanced color eps dl 2.5
filenameoutput="CDA-d_FMS_3SAT_Psol_stepsdec_5_paper.eps"
set output filenameoutput

reset
set size 1,1
set origin 0,0
set ylabel "{/=22 P(sol)}" rotate by 90 offset 1,-1
set xlabel "{/=22 {/Symbol a}}" offset 0,0
set key spacing 2.0 maxrows 10 width 6 at 4.19,1.0
set xrange [3.79:4.21]
set yrange [0:1.05]
set tics  font ",24"
set format y "%g"


ac = 3.86
unset arrow
set arrow from ac,graph 0 to ac,graph 1 nohead dt 2
set label "{/=22 {/Symbol a}_{d}}" at 3.87, 0.05 


p "< awk '$1==512 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 1 lw 2 t '{/=24 N=512',\
"< awk '$1==1024 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 2 lw 2 t '{/=24 N=1024',\
"< awk '$1==2048 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 4 lw 2 t '{/=24 N=2048',\
"< awk '$1==4096 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 6 lw 2 t '{/=24 N=4096',\
"< awk '$1==8192 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 7 lw 2 t '{/=24 N=8192',\
"< awk '$1==16384 && $3==5' CDA_decimation_FMS_all_K_3.txt" u 2:6 w p ps 1.2 lc 8 lw 2 t '{/=24 N=16384',\
"sigmoid_data_N_512_stepsdec_5.txt" u 1:2 w l lc 1 lw 2 notitle,\
"sigmoid_data_N_1024_stepsdec_5.txt" u 1:2 w l lc 2 lw 2 notitle,\
"sigmoid_data_N_2048_stepsdec_5.txt" u 1:2 w l lc 4 lw 2 notitle,\
"sigmoid_data_N_4096_stepsdec_5.txt" u 1:2 w l lc 6 lw 2 notitle,\
"sigmoid_data_N_8192_stepsdec_5.txt" u 1:2 w l lc 7 lw 2 notitle,\
"sigmoid_data_N_16384_stepsdec_5.txt" u 1:2 w l lc 8 lw 2 notitle