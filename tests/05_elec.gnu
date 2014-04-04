set terminal postscript eps color enhanced "Helvetica" 15
set output "05_elec.eps" 

set style line  1 lt 1 lw 1.5 lc rgb "#000000" pt  4 ps 0.75
set style line  2 lt 1 lw 1.5 lc rgb "#FF0000" pt  4 ps 0.75
set style line  3 lt 1 lw 1.5 lc rgb "#0000FF" pt  4 ps 0.75
set style line  4 lt 2 lw 1.5 lc rgb "#000000" pt  8 ps 0.75
set style line  5 lt 2 lw 1.5 lc rgb "#FF0000" pt  8 ps 0.75
set style line  6 lt 2 lw 1.5 lc rgb "#0000FF" pt  8 ps 0.75
set style line  7 lt 4 lw 1.5 lc rgb "#000000" pt  4 ps 0.75
set style line  8 lt 4 lw 1.5 lc rgb "#FF0000" pt  6 ps 0.75
set style line 11 lc rgb "#0000FF" pt 1 ps 1.3 lw 2
set style line 12 lc rgb "#000000" pt 2 ps 1.3 lw 2
set style line 13 lc rgb "#AC8700" pt 12 ps 1 lw 3
set style line 14 lc rgb "#0087AC" pt 10 ps 1 lw 3
set style line 15 lc rgb "#5500FE" pt 14 ps 1 lw 3

set xlabel 'Volume of prim. cell (bohr^3)'
set ylabel 'Free energy of prim. cell (kJ/mol)'
set format y "%.2f"

plot '05_elec.eos' u 3:23 index 1 w linesp ls 1 title 'F_{el}, free electron',\
     '05_elec.eos' u 3:23 index 3 w linesp ls 2 title 'F_{el}, free electron, calc. N(Ef)',\
     '05_elec.eos' u 3:23 index 5 w linesp ls 3 title 'F_{el}, finite T DFT',\
     '05_elec.eos' u 3:(-$27*$2/1000) index 1 w linesp ls 4 title '-T*S_{el}, free electron',\
     '05_elec.eos' u 3:(-$27*$2/1000) index 3 w linesp ls 5 title '-T*S_{el}, free electron, calc. N(Ef)',\
     '05_elec.eos' u 3:(-$27*$2/1000) index 5 w linesp ls 6 title '-T*S_{el}, finite T DFT'
     

!epstopdf 05_elec.eps
!pdfcrop 05_elec.pdf
!mv 05_elec-crop.pdf 05_elec.pdf
!rm 05_elec.eps
