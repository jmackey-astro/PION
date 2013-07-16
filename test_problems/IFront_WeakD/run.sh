#!/bin/bash

rm data/*

../../icgen_serial params_WeakD_n0320.txt silo
../../icgen_serial params_WeakD_n0640.txt silo
../../icgen_serial params_WeakD_n1280.txt silo
../../icgen_serial params_WeakD_n2560.txt silo

../../pion_serial IC_IF_WeakD_MPv7_n0320.silo 5 1 optype=text redirect=data/n0320_ &
../../pion_serial IC_IF_WeakD_MPv7_n0640.silo 5 1 optype=text redirect=data/n0640_ &
../../pion_serial IC_IF_WeakD_MPv7_n1280.silo 5 1 optype=text redirect=data/n1280_ &
../../pion_serial IC_IF_WeakD_MPv7_n2560.silo 5 1 optype=text redirect=data/n2560_ &

wait

exit
