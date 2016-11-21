set terminal postscript eps font 'Helvetica, 24'
set output 'e01d0001.eps'

set xlabel '# samples'
set ylabel 'error bounds'

plot 'bounds1.txt' using 1:2 title 'our bound' with linespoints pt 3, \
	'bounds1.txt' using 1:3 title 'Riondato and Upfal' with linespoints pt 4	
set term epslatex
set output
