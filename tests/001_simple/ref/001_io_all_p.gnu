set terminal postscript eps color enhanced "Helvetica" 14

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

!cp /home/alberto/git/gibbs2/build/tests/001_simple/001_io.eos temp.dat

set xrange [0:               146.5]
set xlabel "p (GPa)"

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03.eps'
set ylabel "V(bohr^3)"
plot \
 'temp.dat' u 1: 3 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 3 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 3 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 3 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 3 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_03.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04.eps'
set ylabel "Estatic(Ha)"
plot \
 'temp.dat' u 1: 4 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 4 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 4 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 4 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 4 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_04.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05.eps'
set ylabel "G(kJ/mol)"
plot \
 'temp.dat' u 1: 5 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 5 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 5 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 5 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 5 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_05.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06.eps'
set ylabel "Gerr(kJ/mol)"
plot \
 'temp.dat' u 1: 6 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 6 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 6 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 6 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 6 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_06.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07.eps'
set ylabel "p_sta(GPa)"
plot \
 'temp.dat' u 1: 7 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 7 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 7 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 7 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 7 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_07.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08.eps'
set ylabel "p_th(GPa)"
plot \
 'temp.dat' u 1: 8 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 8 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 8 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 8 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 8 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_08.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09.eps'
set ylabel "B(GPa)"
plot \
 'temp.dat' u 1: 9 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1: 9 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1: 9 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1: 9 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1: 9 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_09.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10.eps'
set ylabel "U-Esta(kJ/mol)"
plot \
 'temp.dat' u 1:10 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:10 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:10 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:10 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:10 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_10.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11.eps'
set ylabel "Cv(J/molK)"
plot \
 'temp.dat' u 1:11 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:11 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:11 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:11 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:11 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_11.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12.eps'
set ylabel "F-Esta(kJ/mol)"
plot \
 'temp.dat' u 1:12 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:12 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:12 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:12 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:12 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_12.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13.eps'
set ylabel "S(J/molK)"
plot \
 'temp.dat' u 1:13 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:13 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:13 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:13 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:13 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_13.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14.eps'
set ylabel "ThetaD(K)"
plot \
 'temp.dat' u 1:14 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:14 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:14 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:14 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:14 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_14.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15.eps'
set ylabel "gamma"
plot \
 'temp.dat' u 1:15 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:15 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:15 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:15 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:15 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_15.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16.eps'
set ylabel "alpha(10^-5/K)"
plot \
 'temp.dat' u 1:16 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:16 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:16 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:16 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:16 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_16.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17.eps'
set ylabel "dp/dT(GPa/K)"
plot \
 'temp.dat' u 1:17 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:17 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:17 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:17 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:17 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_17.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18.eps'
set ylabel "Bs(GPa)"
plot \
 'temp.dat' u 1:18 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:18 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:18 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:18 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:18 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_18.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19.eps'
set ylabel "Cp(J/molK)"
plot \
 'temp.dat' u 1:19 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:19 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:19 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:19 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:19 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_19.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20.eps'
set ylabel "B_Tp"
plot \
 'temp.dat' u 1:20 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:20 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:20 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:20 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:20 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_20.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21.eps'
set ylabel "B_Tpp(GPa-1)"
plot \
 'temp.dat' u 1:21 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:21 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:21 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:21 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:21 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_21.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22.eps'
set ylabel "Fvib(kJ/mol)"
plot \
 'temp.dat' u 1:22 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:22 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:22 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:22 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:22 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_22.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23.eps'
set ylabel "Fel(kJ/mol)"
plot \
 'temp.dat' u 1:23 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:23 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:23 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:23 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:23 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_23.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24.eps'
set ylabel "Uvib(kJ/mol)"
plot \
 'temp.dat' u 1:24 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:24 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:24 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:24 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:24 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_24.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25.eps'
set ylabel "Uel(kJ/mol)"
plot \
 'temp.dat' u 1:25 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:25 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:25 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:25 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:25 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_25.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26.eps'
set ylabel "Svib(J/molK)"
plot \
 'temp.dat' u 1:26 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:26 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:26 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:26 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:26 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_26.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27.eps'
set ylabel "Sel(J/molK)"
plot \
 'temp.dat' u 1:27 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:27 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:27 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:27 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:27 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_27.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28.eps'
set ylabel "Cv_vib(J/molK)"
plot \
 'temp.dat' u 1:28 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:28 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:28 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:28 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:28 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_28.eps

set output '/home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29.eps'
set ylabel "Cv_el(J/molK)"
plot \
 'temp.dat' u 1:29 index   0 w lines ls  1 title '       mgo,    0.00K'  ,\
 'temp.dat' u 1:29 index  25 w lines ls  2 title '       mgo,  165.08K'  ,\
 'temp.dat' u 1:29 index  50 w lines ls  3 title '       mgo,  330.15K'  ,\
 'temp.dat' u 1:29 index  75 w lines ls  4 title '       mgo,  495.23K'  ,\
 'temp.dat' u 1:29 index  99 w lines ls  5 title '       mgo,  653.70K'
!epstopdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29.eps
!pdfcrop /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29.pdf
!mv /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29-crop.pdf /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29.pdf
!rm /home/alberto/git/gibbs2/build/tests/001_simple/001_io_p_29.eps

!rm temp.dat
