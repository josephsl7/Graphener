#HFE_plot.gp, formation energy referenced to H2 and the metal in a hexagonal monolayer, isolated.
set terminal postscript landscape color "Helvetica" 24 
set output "HFE.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title ""
set nokey
plot "vaspHFE.out" using 2:3 lt 1 pt 6

!ps2pdf HFE.eps && rm HFE.eps

