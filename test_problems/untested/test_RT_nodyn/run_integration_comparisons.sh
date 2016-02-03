#!/bin/bash
#
# File to compare the 2nd order accurate ray-tracing scheme for different timestepping criteria in xdot.
#


##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err050/
mkdir res201106_O2err050/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err050* res201106_O2err050/res/*O2err050*
#mv res201106_O2err050/res/*O2err050* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err050
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err050* res201106_O2err050/res
mv COMP_* res201106_O2err050
#
##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err100/
mkdir res201106_O2err100/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err100* res201106_O2err100/res/*O2err100*
#mv res201106_O2err100/res/*O2err100* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err100
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err100* res201106_O2err100/res
mv COMP_* res201106_O2err100
#
##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err200/
mkdir res201106_O2err200/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err200* res201106_O2err200/res/*O2err200* 
#mv res201106_O2err200/res/*O2err200* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err200
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err200* res201106_O2err200/res
mv COMP_* res201106_O2err200
#
##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err300/
mkdir res201106_O2err300/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err300* res201106_O2err300/res/*O2err300*
#mv res201106_O2err300/res/*O2err300* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err300
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err300* res201106_O2err300/res
mv COMP_* res201106_O2err300
#
##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err400/
mkdir res201106_O2err400/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err400* res201106_O2err400/res/*O2err400*
#mv res201106_O2err400/res/*O2err400* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err400
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err400* res201106_O2err400/res
mv COMP_* res201106_O2err400
#
##
## Run the comparison code for each timestepping criterion in turn.
##
mkdir res201106_O2err500/
mkdir res201106_O2err500/res
rm /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err500* res201106_O2err500/res/*O2err500*
#mv res201106_O2err500/res/*O2err500* /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/
./compare_RT_methods.sh O2err500
mv /vol/aibn129/aibn129_1/jmackey/current_data/code_tests/RT_static_comp/*O2err500* res201106_O2err500/res
mv COMP_* res201106_O2err500
#
##
## Now view the results:
##
#eog  res201106_O2err050/COMP_*.jpeg
#eog  res201106_O2err100/COMP_*.jpeg 
#eog  res201106_O2err200/COMP_*.jpeg
#eog  res201106_O2err300/COMP_*.jpeg
#eog  res201106_O2err400/COMP_*.jpeg
#eog  res201106_O2err500/COMP_*.jpeg

#wait
exit

