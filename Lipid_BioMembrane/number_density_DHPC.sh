#! /bin/bash

#####################################################################
# Calculates the number of DPPC/DHPC lipids within certain radius  ##
#####################################################################

b_0=225000 ## initial tile of extraction of traj file from whole trajectory
e_0=250000
dt=300 ## time step for the trajectory file extraction
i=0
y=8   ##number of beads in lipids
sum=0
count1=0
sq_dev=0
count2=0
z=0.5
#echo 5|trjconv_465_mpi -f traj0.xtc -b $b_0 -e $e_0 -dt $dt -n index.ndx -o traj.xtc
#sleep 1
g_select_465_mpi -sf selection.dat -selrpos whole_res_com  -f traj.xtc -s topol.tpr -on radius.ndx

for ((b=225000 ,e=225300; e<250000; b+=$dt,e+=$dt));
    do
     echo $i|trjconv_465_mpi -f traj.xtc -b $b -e $e -n radius.ndx -o dppc_radius_$i.gro  
     sed -n 2p dppc_radius_$i.gro>>avg_no.dat
     rm dppc_radius_$i.gro 
     i=$(($i+1))
    done

###########################################
#Calculate Averages and Standard Error ####
###########################################

 
 for x in `cat avg_no.dat`;
    do 
     sum=$(($sum + $x))
     count1=$(($count1+1))
    done
     avg=$(echo $sum $count1  | awk '{printf "%.3f",$1 / $2 }')
     no_of_lipids=$(echo $avg $y | awk '{printf "%.3f",$1 / $2}')

 for j in `cat avg_no.dat`;
     do
       dev=$(echo $j $avg | awk '{printf "%.3f", $1 - $2}')
       sq_dev=$(echo $sq_dev $dev  | awk '{printf "%.3f",$1 + ($2 * $2)}')
       count2=$(($count2+1))
     done
       variance=$(echo $sq_dev $count2 | awk '{printf "%.3f",$1 / ($2 - 1)}')
       sd=$(echo $variance | awk '{printf "%.3f", sqrt($1)}')
       se=$(echo $sd $count2 $y | awk '{printf "%.3f", $1 / (sqrt($2) * $3)}')
  
    
     echo "Number of Lipids= "$no_of_lipids" "+/-" "$se" "
     echo " "
     echo "Number of Lipids= "$no_of_lipids" "+/-" "$se" ">>sum.dat
     echo>avg_no.dat

