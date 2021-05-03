#!/usr/local/bin/gnuplot

set terminal postscript eps enhanced background rgb 'white' color size 10cm , 10cm  "Times" 24
set output "Fig2_revised.eps"

set encoding iso_8859_1
set size square
set tics out
set tics font "Times-Roman,24"
set xlabel font "Times-Roman,24"
set ylabel font "Times-Roman,24"
set label font "Times-Roman,24"
set key font "Times-Roman,24"

set style fill solid border lc rgb "black"

set xlabel "^{3}J_{exp.}"

set ylabel "^{3}J_{cal.}"

xmin=2
xmax=11
set xrange [xmin:xmax]

ymin=2
ymax=11
set yrange [ymin:ymax]

f(x)=x

plot "J3_HN_HA_6LYT_100nsave_vs_exp.txt" u 1:2 w p pt 7 ps 2.0 lc rgb "black" notitle, \
     "J3_HN_HA_6LYT_100nswave_vs_exp.txt" u 1:2 w p pt 7 ps 2.0 lc rgb "red" notitle, \
     "J3_HN_HA_6LYT_100nsave_vs_exp_select.txt" u 1:2 w p pt 6 ps 3.0 lc rgb "black" notitle, \
     "J3_HN_HA_6LYT_100nswave_vs_exp_select.txt" u 1:2 w p pt 6 ps 3.0 lc rgb "black" notitle, \
     f(x) w l lt 1 lc rgb "black" notitle

quit
