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

!sed '/^ *$/d' 001_fits.eos | awk '$1+0==0' | awk '/# Phase/{print ""; print ""} {print}'> temp.dat

set xrange [0:              3000.0]
set xlabel "T (K)"

set output '001_fits_t_03.eps'
set ylabel "V(bohr^3)"
plot \
     'temp.dat' u 2: 3 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 3 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 3 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 3 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 3 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_03.eps
!pdfcrop 001_fits_t_03.pdf
!mv 001_fits_t_03-crop.pdf 001_fits_t_03.pdf
!rm 001_fits_t_03.eps

set output '001_fits_t_04.eps'
set ylabel "Estatic(Ha)"
plot \
     'temp.dat' u 2: 4 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 4 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 4 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 4 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 4 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_04.eps
!pdfcrop 001_fits_t_04.pdf
!mv 001_fits_t_04-crop.pdf 001_fits_t_04.pdf
!rm 001_fits_t_04.eps

set output '001_fits_t_05.eps'
set ylabel "G(kJ/mol)"
plot \
     'temp.dat' u 2: 5 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 5 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 5 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 5 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 5 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_05.eps
!pdfcrop 001_fits_t_05.pdf
!mv 001_fits_t_05-crop.pdf 001_fits_t_05.pdf
!rm 001_fits_t_05.eps

set output '001_fits_t_06.eps'
set ylabel "Gerr(kJ/mol)"
plot \
     'temp.dat' u 2: 6 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 6 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 6 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 6 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 6 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_06.eps
!pdfcrop 001_fits_t_06.pdf
!mv 001_fits_t_06-crop.pdf 001_fits_t_06.pdf
!rm 001_fits_t_06.eps

set output '001_fits_t_07.eps'
set ylabel "p_sta(GPa)"
plot \
     'temp.dat' u 2: 7 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 7 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 7 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 7 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 7 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_07.eps
!pdfcrop 001_fits_t_07.pdf
!mv 001_fits_t_07-crop.pdf 001_fits_t_07.pdf
!rm 001_fits_t_07.eps

set output '001_fits_t_08.eps'
set ylabel "p_th(GPa)"
plot \
     'temp.dat' u 2: 8 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 8 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 8 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 8 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 8 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_08.eps
!pdfcrop 001_fits_t_08.pdf
!mv 001_fits_t_08-crop.pdf 001_fits_t_08.pdf
!rm 001_fits_t_08.eps

set output '001_fits_t_09.eps'
set ylabel "B(GPa)"
plot \
     'temp.dat' u 2: 9 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2: 9 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2: 9 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2: 9 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2: 9 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_09.eps
!pdfcrop 001_fits_t_09.pdf
!mv 001_fits_t_09-crop.pdf 001_fits_t_09.pdf
!rm 001_fits_t_09.eps

set output '001_fits_t_10.eps'
set ylabel "U-Esta(kJ/mol)"
plot \
     'temp.dat' u 2:10 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:10 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:10 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:10 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:10 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_10.eps
!pdfcrop 001_fits_t_10.pdf
!mv 001_fits_t_10-crop.pdf 001_fits_t_10.pdf
!rm 001_fits_t_10.eps

set output '001_fits_t_11.eps'
set ylabel "Cv(J/molK)"
plot \
     'temp.dat' u 2:11 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:11 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:11 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:11 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:11 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_11.eps
!pdfcrop 001_fits_t_11.pdf
!mv 001_fits_t_11-crop.pdf 001_fits_t_11.pdf
!rm 001_fits_t_11.eps

set output '001_fits_t_12.eps'
set ylabel "F-Esta(kJ/mol)"
plot \
     'temp.dat' u 2:12 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:12 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:12 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:12 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:12 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_12.eps
!pdfcrop 001_fits_t_12.pdf
!mv 001_fits_t_12-crop.pdf 001_fits_t_12.pdf
!rm 001_fits_t_12.eps

set output '001_fits_t_13.eps'
set ylabel "S(J/molK)"
plot \
     'temp.dat' u 2:13 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:13 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:13 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:13 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:13 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_13.eps
!pdfcrop 001_fits_t_13.pdf
!mv 001_fits_t_13-crop.pdf 001_fits_t_13.pdf
!rm 001_fits_t_13.eps

set output '001_fits_t_14.eps'
set ylabel "ThetaD(K)"
plot \
     'temp.dat' u 2:14 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:14 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:14 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:14 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:14 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_14.eps
!pdfcrop 001_fits_t_14.pdf
!mv 001_fits_t_14-crop.pdf 001_fits_t_14.pdf
!rm 001_fits_t_14.eps

set output '001_fits_t_15.eps'
set ylabel "gamma"
plot \
     'temp.dat' u 2:15 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:15 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:15 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:15 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:15 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_15.eps
!pdfcrop 001_fits_t_15.pdf
!mv 001_fits_t_15-crop.pdf 001_fits_t_15.pdf
!rm 001_fits_t_15.eps

set output '001_fits_t_16.eps'
set ylabel "alpha(10^-5/K)"
plot \
     'temp.dat' u 2:16 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:16 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:16 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:16 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:16 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_16.eps
!pdfcrop 001_fits_t_16.pdf
!mv 001_fits_t_16-crop.pdf 001_fits_t_16.pdf
!rm 001_fits_t_16.eps

set output '001_fits_t_17.eps'
set ylabel "dp/dT(GPa/K)"
plot \
     'temp.dat' u 2:17 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:17 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:17 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:17 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:17 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_17.eps
!pdfcrop 001_fits_t_17.pdf
!mv 001_fits_t_17-crop.pdf 001_fits_t_17.pdf
!rm 001_fits_t_17.eps

set output '001_fits_t_18.eps'
set ylabel "Bs(GPa)"
plot \
     'temp.dat' u 2:18 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:18 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:18 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:18 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:18 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_18.eps
!pdfcrop 001_fits_t_18.pdf
!mv 001_fits_t_18-crop.pdf 001_fits_t_18.pdf
!rm 001_fits_t_18.eps

set output '001_fits_t_19.eps'
set ylabel "Cp(J/molK)"
plot \
     'temp.dat' u 2:19 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:19 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:19 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:19 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:19 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_19.eps
!pdfcrop 001_fits_t_19.pdf
!mv 001_fits_t_19-crop.pdf 001_fits_t_19.pdf
!rm 001_fits_t_19.eps

set output '001_fits_t_20.eps'
set ylabel "B_Tp"
plot \
     'temp.dat' u 2:20 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:20 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:20 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:20 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:20 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_20.eps
!pdfcrop 001_fits_t_20.pdf
!mv 001_fits_t_20-crop.pdf 001_fits_t_20.pdf
!rm 001_fits_t_20.eps

set output '001_fits_t_21.eps'
set ylabel "B_Tpp(GPa-1)"
plot \
     'temp.dat' u 2:21 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:21 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:21 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:21 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:21 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_21.eps
!pdfcrop 001_fits_t_21.pdf
!mv 001_fits_t_21-crop.pdf 001_fits_t_21.pdf
!rm 001_fits_t_21.eps

set output '001_fits_t_22.eps'
set ylabel "Fvib(kJ/mol)"
plot \
     'temp.dat' u 2:22 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:22 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:22 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:22 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:22 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_22.eps
!pdfcrop 001_fits_t_22.pdf
!mv 001_fits_t_22-crop.pdf 001_fits_t_22.pdf
!rm 001_fits_t_22.eps

set output '001_fits_t_23.eps'
set ylabel "Fel(kJ/mol)"
plot \
     'temp.dat' u 2:23 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:23 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:23 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:23 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:23 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_23.eps
!pdfcrop 001_fits_t_23.pdf
!mv 001_fits_t_23-crop.pdf 001_fits_t_23.pdf
!rm 001_fits_t_23.eps

set output '001_fits_t_24.eps'
set ylabel "Uvib(kJ/mol)"
plot \
     'temp.dat' u 2:24 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:24 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:24 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:24 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:24 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_24.eps
!pdfcrop 001_fits_t_24.pdf
!mv 001_fits_t_24-crop.pdf 001_fits_t_24.pdf
!rm 001_fits_t_24.eps

set output '001_fits_t_25.eps'
set ylabel "Uel(kJ/mol)"
plot \
     'temp.dat' u 2:25 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:25 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:25 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:25 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:25 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_25.eps
!pdfcrop 001_fits_t_25.pdf
!mv 001_fits_t_25-crop.pdf 001_fits_t_25.pdf
!rm 001_fits_t_25.eps

set output '001_fits_t_26.eps'
set ylabel "Svib(J/molK)"
plot \
     'temp.dat' u 2:26 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:26 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:26 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:26 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:26 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_26.eps
!pdfcrop 001_fits_t_26.pdf
!mv 001_fits_t_26-crop.pdf 001_fits_t_26.pdf
!rm 001_fits_t_26.eps

set output '001_fits_t_27.eps'
set ylabel "Sel(J/molK)"
plot \
     'temp.dat' u 2:27 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:27 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:27 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:27 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:27 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_27.eps
!pdfcrop 001_fits_t_27.pdf
!mv 001_fits_t_27-crop.pdf 001_fits_t_27.pdf
!rm 001_fits_t_27.eps

set output '001_fits_t_28.eps'
set ylabel "Cv_vib(J/molK)"
plot \
     'temp.dat' u 2:28 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:28 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:28 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:28 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:28 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_28.eps
!pdfcrop 001_fits_t_28.pdf
!mv 001_fits_t_28-crop.pdf 001_fits_t_28.pdf
!rm 001_fits_t_28.eps

set output '001_fits_t_29.eps'
set ylabel "Cv_el(J/molK)"
plot \
     'temp.dat' u 2:29 index   0 w lines title '       bm3'  ,\
     'temp.dat' u 2:29 index   1 w lines title '       bm4'  ,\
     'temp.dat' u 2:29 index   2 w lines title '       pt3'  ,\
     'temp.dat' u 2:29 index   3 w lines title '       pt4'  ,\
     'temp.dat' u 2:29 index   4 w lines title '     vinet'
!epstopdf 001_fits_t_29.eps
!pdfcrop 001_fits_t_29.pdf
!mv 001_fits_t_29-crop.pdf 001_fits_t_29.pdf
!rm 001_fits_t_29.eps

!rm temp.dat
