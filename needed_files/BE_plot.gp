# BE_plot.gp, Binding energy vs. graphene
set terminal postscript landscape color "Helvetica" 24 
set output "BE.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title ""
set nokey
plot "vaspBE.out" using 2:3 lt 1 pt 6

!ps2pdf BE.eps && rm BE.eps

