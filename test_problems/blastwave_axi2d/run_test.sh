#!/bin/bash

mpirun -np 16 ../../pion-ng-bc4 BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa2bc4 opfreq=1024 redirect=log_np04oa2bc4 ooa=2 &
mpirun -np 16 ../../pion-ng-bc6 BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa2bc6 opfreq=1024 redirect=log_np04oa2bc6 ooa=2 &
mpirun -np 16 ../../pion-ng-cf2 BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa2cf2 opfreq=1024 redirect=log_np04oa2cf2 ooa=2 &
wait
exit

#mpirun -np 32 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np32oa2 opfreq=1024 redirect=log_np32oa2 ooa=2 finishtime=5.0e10 &
#mpirun -np 16 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np16oa2 opfreq=1024 redirect=log_np16oa2 ooa=2 finishtime=5.0e10 &
mpirun -np  8 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np08oa2bc4 opfreq=1024 redirect=log_np08oa2bc4 ooa=2 finishtime=5.0e10 
mpirun -np  4 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa2bc4 opfreq=1024 redirect=log_np04oa2bc4 ooa=2 finishtime=5.0e10 &
mpirun -np  2 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np02oa2bc4 opfreq=1024 redirect=log_np02oa2bc4 ooa=2 finishtime=5.0e10 &
mpirun -np  1 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np01oa2bc4 opfreq=1024 redirect=log_np01oa2bc4 ooa=2 finishtime=5.0e10 &
wait

grep TOTALS log_np*.txt
grep "New time: 5.0000000e+10" log_np*
exit

#mpirun -np 32 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np32oa2 opfreq=1024 redirect=log_np32oa2 ooa=2 finishtime=5.0e10 &
#mpirun -np 16 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np16oa2 opfreq=1024 redirect=log_np16oa2 ooa=2 finishtime=5.0e10 &
mpirun -np  8 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np08oa2cf opfreq=1024 redirect=log_np08oa2cf ooa=2 finishtime=5.0e10 
mpirun -np  4 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa2cf opfreq=1024 redirect=log_np04oa2cf ooa=2 finishtime=5.0e10 &
mpirun -np  2 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np02oa2cf opfreq=1024 redirect=log_np02oa2cf ooa=2 finishtime=5.0e10 &
mpirun -np  1 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np01oa2cf opfreq=1024 redirect=log_np01oa2cf ooa=2 finishtime=5.0e10 &
wait

grep TOTALS log_np*.txt
grep "New time: 5.0000000e+10" log_np*
exit

mpirun -np 32 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np32oa1 opfreq=1024 redirect=log_np32oa1 ooa=1 finishtime=2.0e11 &
mpirun -np 16 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np16oa1 opfreq=1024 redirect=log_np16oa1 ooa=1 finishtime=2.0e11 &
mpirun -np  8 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np08oa1 opfreq=1024 redirect=log_np08oa1 ooa=1 finishtime=2.0e11 &
mpirun -np  4 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04oa1 opfreq=1024 redirect=log_np04oa1 ooa=1 finishtime=2.0e11 &
mpirun -np  2 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np02oa1 opfreq=1024 redirect=log_np02oa1 ooa=1 finishtime=2.0e11 &
mpirun -np  1 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np01oa1 opfreq=1024 redirect=log_np01oa1 ooa=1 finishtime=2.0e11 &
wait

grep TOTALS log_np*.txt
exit

mpirun -np 32 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np32 opfreq=1024 redirect=log_np32 &
mpirun -np 16 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np16 opfreq=1024 redirect=log_np16 &
mpirun -np  8 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np08 opfreq=1024 redirect=log_np08 &
mpirun -np  4 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np04 opfreq=1024 redirect=log_np04 &
mpirun -np  2 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np02 opfreq=1024 redirect=log_np02 &
mpirun -np  1 ../../pion-ng BW_NGaxi2D_NR128_level00_0000.00000000.silo outfile=NGaxi2D_NR128_np01 opfreq=1024 redirect=log_np01 &
wait

grep TOTALS log_np*.txt
exit

