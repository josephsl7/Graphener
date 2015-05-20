# BE_plot.gp, Binding energy vs. graphene
set terminal postscript landscape color "Helvetica" 24 
set output "BE.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title ""
 
plot "uncleBE.out" using 2:4 lt 1 lw 3 pt 2 lc rgb 'black' title "Fit/Pred" , "vaspBE.out" using 2:3 lt 1 lw 3  pt 65 lc rgb 'red' title "VASP"

!ps2pdf BE.eps && rm BE.eps

