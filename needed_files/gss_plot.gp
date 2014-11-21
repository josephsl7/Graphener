# gss_plot.gp
set terminal postscript landscape color "Helvetica" 24 
set output "gss.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title "Ground State Search"
set nokey
plot "gss.out" using 3:8 lt 3, "vaspFE.out" using 2:3 lt 1 pt 6, "gssFailedVasp.out" using 2:3 lt 2

!ps2pdf gss.eps && rm gss.eps

