#! /bin/bash
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=0 and rdist<6/resname DPPC and rdist>=6 and rdist<8/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=6 and rdist<8/resname DPPC and rdist>=8 and rdist<10/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=8 and rdist<10/resname DPPC and rdist>=10 and rdist<12/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=10 and rdist<12/resname DPPC and rdist>=12 and rdist<13/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=12 and rdist<13/resname DPPC and rdist>=13 and rdist<14/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=13 and rdist<14/resname DPPC and rdist>=14 and rdist<15/g' selection.dat
./number_density_DPPC.sh
 sed -i 's/resname DPPC and rdist>=14 and rdist<15/resname DPPC and rdist>=15/g' selection.dat
./number_density_DPPC.sh

########################################################################
 sed -i 's/resname DPPC and rdist>=15/resname DHPC and rdist>=0 and rdist<6/g' selection.dat
########################################################################
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=0 and rdist<6/resname DHPC and rdist>=6 and rdist<8/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=6 and rdist<8/resname DHPC and rdist>=8 and rdist<10/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=8 and rdist<10/resname DHPC and rdist>=10 and rdist<12/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=10 and rdist<12/resname DHPC and rdist>=12 and rdist<13/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=12 and rdist<13/resname DHPC and rdist>=13 and rdist<14/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=13 and rdist<14/resname DHPC and rdist>=14 and rdist<15/g' selection.dat
./number_density_DHPC.sh
 sed -i 's/resname DHPC and rdist>=14 and rdist<15/resname DHPC and rdist>=15/g' selection.dat
./number_density_DHPC.sh
##############################################################################
 sed -i 's/resname DHPC and rdist>=15/resname DPPC and rdist>=0 and rdist<6/g' selection.dat
##############################################################################

