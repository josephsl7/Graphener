# gss_plot.gp
set terminal postscript landscape color "Helvetica" 24 
set output "gss.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title "Ground State Search"

plot "gss.out" using 2:8 lt 1 lw 3 pt 2 lc rgb 'black' title "Fit/Pred", "vaspFE.out" using 2:3 lt 1 lw 3  pt 65 lc rgb 'red' title "VASP", "gssFailedVasp.out" using 2:3 lt 1 lw 3  pt 12 title "Failed"

!ps2pdf gss.eps && rm gss.eps

