#!/bin/bash


for eta in 1.0
do

echo "#!/bin/bash" > "scr_Weigt_KSAT_eta_"$eta".sh"
t=10

for alpha in 2.7 2.71 2.72
do

echo "python Weigt_KSAT_cluster.py "$eta" "$alpha" "$t"  > out_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt 2> err_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt" >> "scr_Weigt_KSAT_eta_"$eta".sh"

echo "done eta="$eta"  alpha="$alpha""

done

chmod +x ./"scr_Weigt_KSAT_eta_"$eta".sh"
nohup ./"scr_Weigt_KSAT_eta_"$eta".sh" > "out_Weigt_KSAT_eta_"$eta".txt" 2> "err_Weigt_KSAT_eta_"$eta".txt" & 

done


###############################################################################################################


for eta in 0.9
do

echo "#!/bin/bash" > "scr_Weigt_KSAT_eta_"$eta".sh"
t=10

for alpha in 2.72 2.75 2.78
do

echo "python Weigt_KSAT_cluster.py "$eta" "$alpha" "$t"  > out_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt 2> err_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt" >> "scr_Weigt_KSAT_eta_"$eta".sh"

echo "done eta="$eta"  alpha="$alpha""

done

chmod +x ./"scr_Weigt_KSAT_eta_"$eta".sh"
nohup ./"scr_Weigt_KSAT_eta_"$eta".sh" > "out_Weigt_KSAT_eta_"$eta".txt" 2> "err_Weigt_KSAT_eta_"$eta".txt" & 

done


###############################################################################################################


for eta in 0.8
do

echo "#!/bin/bash" > "scr_Weigt_KSAT_eta_"$eta".sh"
t=10

for alpha in 2.75 2.78 2.8
do

echo "python Weigt_KSAT_cluster.py "$eta" "$alpha" "$t"  > out_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt 2> err_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt" >> "scr_Weigt_KSAT_eta_"$eta".sh"

echo "done eta="$eta"  alpha="$alpha""

done

chmod +x ./"scr_Weigt_KSAT_eta_"$eta".sh"
nohup ./"scr_Weigt_KSAT_eta_"$eta".sh" > "out_Weigt_KSAT_eta_"$eta".txt" 2> "err_Weigt_KSAT_eta_"$eta".txt" & 

done


###############################################################################################################


for eta in 0.7
do

echo "#!/bin/bash" > "scr_Weigt_KSAT_eta_"$eta".sh"
t=10

for alpha in 2.78 2.8 2.85
do

echo "python Weigt_KSAT_cluster.py "$eta" "$alpha" "$t"  > out_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt 2> err_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt" >> "scr_Weigt_KSAT_eta_"$eta".sh"

echo "done eta="$eta"  alpha="$alpha""

done

chmod +x ./"scr_Weigt_KSAT_eta_"$eta".sh"
nohup ./"scr_Weigt_KSAT_eta_"$eta".sh" > "out_Weigt_KSAT_eta_"$eta".txt" 2> "err_Weigt_KSAT_eta_"$eta".txt" & 

done


###############################################################################################################


for eta in 0.6
do

echo "#!/bin/bash" > "scr_Weigt_KSAT_eta_"$eta".sh"
t=10

for alpha in 2.8 2.85 2.9
do

echo "python Weigt_KSAT_cluster.py "$eta" "$alpha" "$t"  > out_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt 2> err_Weigt_KSAT_eta_"$eta"_alpha_"$alpha".txt" >> "scr_Weigt_KSAT_eta_"$eta".sh"

echo "done eta="$eta"  alpha="$alpha""

done

chmod +x ./"scr_Weigt_KSAT_eta_"$eta".sh"
nohup ./"scr_Weigt_KSAT_eta_"$eta".sh" > "out_Weigt_KSAT_eta_"$eta".txt" 2> "err_Weigt_KSAT_eta_"$eta".txt" & 

done