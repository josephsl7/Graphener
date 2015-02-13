# gss_plot.gp
set terminal postscript landscape color "Helvetica" 24 
set output "gss.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title "Ground State Search"

plot "gss.out" using 3:10 lt 3 lw 2 title "Fit" , "vaspFE.out" using 2:3 lt 1 lw 2  pt 12 title "VASP"  #, "gssFailedVasp.out" using 3:2 lt rgb "violet"

!ps2pdf gss.eps && rm gss.eps

