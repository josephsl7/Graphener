# BE_plot.gp, Binding energy vs. graphene
set terminal postscript landscape color "Helvetica" 24 
set output "BE.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title ""
 
plot "vaspBE.out" using 2:3 lt 1 lw 2  pt 12 title "VASP"

!ps2pdf BE.eps && rm BE.eps

