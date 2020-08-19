set terminal postscript eps color enhanced "Helvetica" 14
set output "/home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans.eps"

set style line  1 lt 1 lc rgb "#000000" pt  4 ps 0.75
set style line  2 lt 1 lc rgb "#0000FF" pt  6 ps 0.75
set style line  3 lt 1 lc rgb "#008B00" pt  8 ps 0.75
set style line  4 lt 1 lc rgb "#FF0000" pt 10 ps 0.75
set style line  5 lt 1 lc rgb "#FF00FF" pt 12 ps 0.75
set style line  6 lt 1 lc rgb "#643700" pt 14 ps 0.75
set style line  7 lt 1 lc rgb "#787878" pt  3 ps 0.75
set style line  8 lt 1 lc rgb "#00FFFF" pt  5 ps 0.75
set style line  9 lt 1 lc rgb "#FFB45A" pt  7 ps 0.75
set style line 10 lt 1 lc rgb "#7D26CD" pt  9 ps 0.75
set style line 11 lt 1 lc rgb "#CD9B9B" pt 11 ps 0.75
set style line 12 lt 1 lc rgb "#CD6D0C" pt 13 ps 0.75
set style line 13 lt 1 lc rgb "#00B98B" pt  1 ps 0.75
set style increment user 

set ylabel "p (GPa)"
set xlabel "T (K)"
set yrange [0:            600.0000]
set xrange [0:           1000.0000]
unset key
plot '/home/alberto/git/gibbs2/build/tests/006_phases/001_phases.ptrans' u 1:4 index 0 w linespoints ls 1   
!epstopdf /home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans.pdf
!mv /home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans-crop.pdf /home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans.pdf
!rm /home/alberto/git/gibbs2/build/tests/006_phases/001_phases_ptrans.eps
