#!/bin/bash



rm data/IF_WeakD_Cartesian_MPv7*
../../icgen_serial params_WeakD_Cartesian_n0040.txt silo
../../icgen_serial params_WeakD_Cartesian_n0080.txt silo
../../icgen_serial params_WeakD_Cartesian_n0160.txt silo
../../icgen_serial params_WeakD_Cartesian_n0320.txt silo
../../icgen_serial params_WeakD_Cartesian_n0640.txt silo

../../pion_serial IC_IF_WeakD_Cartesian_MPv7_n0040.silo 5 1 optype=text redirect=data/Crt_n0040_ &
../../pion_serial IC_IF_WeakD_Cartesian_MPv7_n0080.silo 5 1 optype=text redirect=data/Crt_n0080_ &
../../pion_serial IC_IF_WeakD_Cartesian_MPv7_n0160.silo 5 1 optype=text redirect=data/Crt_n0160_ &
../../pion_serial IC_IF_WeakD_Cartesian_MPv7_n0320.silo 5 1 optype=text redirect=data/Crt_n0320_ &
../../pion_serial IC_IF_WeakD_Cartesian_MPv7_n0640.silo 5 1 optype=text redirect=data/Crt_n0640_ &
wait
exit

rm data/IF_WeakD_MPv7*

../../icgen_serial params_WeakD_n0010.txt silo
../../icgen_serial params_WeakD_n0020.txt silo
../../icgen_serial params_WeakD_n0040.txt silo
../../icgen_serial params_WeakD_n0080.txt silo
../../icgen_serial params_WeakD_n0160.txt silo
../../icgen_serial params_WeakD_n0320.txt silo
../../icgen_serial params_WeakD_n0640.txt silo
../../icgen_serial params_WeakD_n1280.txt silo
../../icgen_serial params_WeakD_n2560.txt silo

../../pion_serial IC_IF_WeakD_MPv7_n0010.silo 5 1 optype=text redirect=data/n0010_ &
../../pion_serial IC_IF_WeakD_MPv7_n0020.silo 5 1 optype=text redirect=data/n0020_ &
../../pion_serial IC_IF_WeakD_MPv7_n0040.silo 5 1 optype=text redirect=data/n0040_ &
../../pion_serial IC_IF_WeakD_MPv7_n0080.silo 5 1 optype=text redirect=data/n0080_ &
../../pion_serial IC_IF_WeakD_MPv7_n0160.silo 5 1 optype=text redirect=data/n0160_ &
../../pion_serial IC_IF_WeakD_MPv7_n0320.silo 5 1 optype=text redirect=data/n0320_ &
../../pion_serial IC_IF_WeakD_MPv7_n0640.silo 5 1 optype=text redirect=data/n0640_ &
../../pion_serial IC_IF_WeakD_MPv7_n1280.silo 5 1 optype=text redirect=data/n1280_ &
../../pion_serial IC_IF_WeakD_MPv7_n2560.silo 5 1 optype=text redirect=data/n2560_ &

wait
exit
