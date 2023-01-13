#! /bin/bash
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
 sed -i 's/resname DPPC and rdist>=0 and rdist<6/resname DPPC and rdist>=6 and rdist<8/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
 sed -i 's/resname DPPC and rdist>=6 and rdist<8/resname DPPC and rdist>=8 and rdist<10/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
 sed -i 's/resname DPPC and rdist>=8 and rdist<10/resname DPPC and rdist>=10 and rdist<12/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
 sed -i 's/resname DPPC and rdist>=10 and rdist<12/resname DPPC and rdist>=12 and rdist<14/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
 sed -i 's/resname DPPC and rdist>=12 and rdist<14/resname DPPC and rdist>=14/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DPPC
########################################################################
 sed -i 's/resname DPPC and rdist>=14/resname DHPC and rdist>=0 and rdist<6/g' selection.dat
########################################################################
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
 sed -i 's/resname DHPC and rdist>=0 and rdist<6/resname DHPC and rdist>=6 and rdist<8/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
 sed -i 's/resname DHPC and rdist>=6 and rdist<8/resname DHPC and rdist>=8 and rdist<10/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
 sed -i 's/resname DHPC and rdist>=8 and rdist<10/resname DHPC and rdist>=10 and rdist<12/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
 sed -i 's/resname DHPC and rdist>=10 and rdist<12/resname DHPC and rdist>=12 and rdist<14/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
 sed -i 's/resname DHPC and rdist>=12 and rdist<14/resname DHPC and rdist>=14/g' selection.dat
./do-order-multi-radius.py traj.xtc 475000 500000 400 0 0 1 DHPC
##############################################################################
 sed -i 's/resname DHPC and rdist>=14/resname DPPC and rdist>=0 and rdist<6/g' selection.dat
##############################################################################

