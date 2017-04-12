px=4000
vanish = 0.0001
set style arrow 1 nohead lw 1 lc palette
set palette defined ( 0 "#ffc0c0", vanish "#ffc0c0", vanish "#0000e0", 1 "#00ff00" )

#------- compute ranges -------
set output '/dev/null'
plot 'cur.dat' u 1:2:($3-$1):($4-$2):5 title '' with vectors arrowstyle 1
minv = vanish/GPVAL_CB_MAX
set xrange [GPVAL_X_MIN:GPVAL_X_MAX]
set yrange [GPVAL_Y_MIN:GPVAL_Y_MAX]
set size ratio -1
set autoscale fix

#------- fit 'volt.dat' onto a grid -------
set table 'tmp.txt'
set dgrid 50, 50, 2
splot 'volt.dat' u 1:2:3
unset table

#------- select terminal -------
#-set term png size px,px/(GPVAL_X_MAX-GPVAL_X_MIN)*(GPVAL_Y_MAX-GPVAL_Y_MIN) fontscale px/1000 font "Times,24" lw 5 ; filename='cur.png'
#-set term postscript ; filename='cur.eps'
#-set term pdfcairo size px,px/(GPVAL_X_MAX-GPVAL_X_MIN)*(GPVAL_Y_MAX-GPVAL_Y_MIN) fontscale px/1000 font "Times,24" lw 5 ; filename='cur.pdf'
set term pngcairo size px,px/(GPVAL_X_MAX-GPVAL_X_MIN)*(GPVAL_Y_MAX-GPVAL_Y_MIN) fontscale px/1000 font "Times,24" lw 5 ; filename='cur.png'

#------- plot -------

set output filename

set multiplot

unset cbtics
set cbrange [0:]
#plot "tmp.txt" u 1:2:($3 <= minv ? minv : $3) title '' w image
plot 'volt.dat' u 1:2:($3 <= minv ? minv : $3) title '' with points lw 2 lc palette

set cbtics
plot 'cur.dat' u 1:2:($3-$1):($4-$2):5 title '' with vectors arrowstyle 1

unset multiplot
