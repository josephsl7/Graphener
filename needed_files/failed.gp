# gss_plot.gp
set terminal postscript landscape color "Helvetica" 24 
set output "gss.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title "Ground State Search"
set nokey
plot "gssFailedVasp.out" using 2:3 lt 3

!ps2pdf gssFailedVasp.eps && rm gssFailedVasp.eps

