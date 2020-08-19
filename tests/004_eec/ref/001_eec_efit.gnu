set terminal postscript eps color enhanced "Helvetica" 14
set output "/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.eps"

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

set xrange [             40.0000:             90.0000]
set yrange [            -11.4436:            -11.0535]
set ylabel "Energy (Ha)"
set xlabel "Volume (bohr^3)"
plot '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  0 w points ls  1 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  0 w lines ls  1 title 'lda:uncorr' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  1 w points ls  2 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  1 w lines ls  2 title 'pbe:uncorr' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  2 w points ls  3 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  2 w lines ls  3 title 'lda:pshift' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  3 w points ls  4 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  3 w lines ls  4 title 'pbe:pshift' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  4 w points ls  5 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  4 w lines ls  5 title ' lda:apbaf' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  5 w points ls  6 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  5 w lines ls  6 title ' pbe:apbaf' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  6 w points ls  7 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  6 w lines ls  7 title 'lda:bpscal' ,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec.efit' u 1:2 index  7 w points ls  8 notitle,\
     '/home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.aux' u 1:2 index  7 w lines ls  8 title 'pbe:bpscal'   
!epstopdf /home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.pdf
!mv /home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit-crop.pdf /home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.pdf
!rm /home/alberto/git/gibbs2/build/tests/004_eec/001_eec_efit.eps
