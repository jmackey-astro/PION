#!/bin/bash


simfile=$1
cat << EOF  > gnu.plt
set size 1.0,0.75
#set xrange [2.4e16:3.4e16]
set yrange [-0.1:1.5]
#set yrange [-0.1:1.1]
#set yrange [0.9e-23:1.5e-21]
#set log y
EOF

LIST=`ls ${simfile}*`
for i in $LIST
do
    cat << EOF  >> gnu.plt
#plot '${i}' u 1:2 w l lt -1 
set title "File: ${i}"
plot '${i}' u 1:(log10(\$2)+23.0) w l title "log10(density*1e23)", \
 '${i}' u 1:10 w l title "x(H+)", \
 '${i}' u 1:(1-\$10) w l title "x(H0)", \
 '${i}' u 1:(\$11/1.0e5) w l lt -1 title "T(10^5K)"
#plot '${i}' u 1:12 w l title "He(1+)", \
# '${i}' u 1:13 w l title "He(2+)", \
# '${i}' u 1:(1-\$12-\$13) w l title "He(0)", \
# '${i}' u 1:(\$14/1.0e5) w l lt -1 title "T(10^5K)"
pause -1
EOF
done

gnuplot gnu.plt
