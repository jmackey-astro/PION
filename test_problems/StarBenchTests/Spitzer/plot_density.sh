
DDIR=/vol/klaipeda3/scratch/jmackey/Spitzer_ascii

n0080=`ls ${DDIR}/S1D_n0080.*.txt | tail -n41 | head -n1`
n0160=`ls ${DDIR}/S1D_n0160.*.txt | tail -n41 | head -n1`
n0320=`ls ${DDIR}/S1D_n0320.*.txt | tail -n41 | head -n1`
n0640=`ls ${DDIR}/S1D_n0640.*.txt | tail -n41 | head -n1`
n1280=`ls ${DDIR}/S1D_n1280.*.txt | tail -n41 | head -n1`
n2560=`ls ${DDIR}/S1D_n2560.*.txt | tail -n41 | head -n1`



set log y
var=2
plot "${n0080}" u (\$1/3.086e18):var w lp, \
     "${n0160}" u (\$1/3.086e18):var w lp, \
     "${n0320}" u (\$1/3.086e18):var w lp, \
     "${n0640}" u (\$1/3.086e18):var w l, \
     "${n1280}" u (\$1/3.086e18):var w l, \
     "${n2560}" u (\$1/3.086e18):var w l



