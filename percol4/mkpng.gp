px=4000 ; set term png size px,px/4*3 fontscale px/1000 font "Times,24" lw 5; set output 'cur.png'

set style arrow 1 nohead lw 2 lc palette

set palette rgbformulae 7,5,15
plot 'cur.dat' u 1:2:($3-$1):($4-$2):5 title '' with vectors arrowstyle 1
 