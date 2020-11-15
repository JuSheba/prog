set terminal png size 1000,400
set output 'uniform.png'
plot 'result.dat' using 1:2 with lines
