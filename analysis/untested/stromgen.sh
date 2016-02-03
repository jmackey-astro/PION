
file=../results/pos_pec3drecomb
counter=100
step=10

#for i in $( ls $file* ); 
#do
src=$file.$counter.txt
cat << EOF  > gnu.plt
a=1.e15
f(x)=1-exp(-(a/x)**15)
set yrange [-0.05:1.05]
set xrange [1.e14:3.e16]
fit f(x) '$src' u 2:3:4 via a
plot f(x) title "fit", '$src' u 2:3 title "Sim output"
pause -1
EOF
gnuplot -noraise gnu.plt
rm gnu.plt
echo The counter is $counter
let counter=counter+$step
#done
grep "+/-" fit.log > tmp.txt
sed '/a               = /s///g' tmp.txt  > tmp2.txt
