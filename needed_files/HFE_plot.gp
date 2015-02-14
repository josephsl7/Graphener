#HFE_plot.gp, formation energy referenced to H2 and the metal in a hexagonal monolayer, isolated.
set terminal postscript landscape color "Helvetica" 24 
set output "HFE.eps"
set xlabel "Metal concentration x"
set ylabel "Energy (eV)"
set title ""

plot "uncleHFE.out" using 2:5 lt 1 lw 3 pt 2 lc rgb 'black' title "Fit/Pred" , "vaspHFE.out" using 2:3 lt 1 lw 3  pt 65 lc rgb 'red' title "VASP"

!ps2pdf HFE.eps && rm HFE.eps

