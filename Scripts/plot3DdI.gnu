fileEl='ExVxdIel.tot'
fileTot='ExVxdItot.tot'
fileSy='ExVxdIsy.tot'
fileAsy='ExVxdIasy.tot'
labelCB='dI/dV (A/Ry)'
set xrange [-0.2:0.2]
set yrange [-0.2:0.2]
set cbrange [-5e-04:5e-04]

set terminal jpeg large size 2200, 2000 enhanced font "Arial,20"

#set xtics font "Arial,20" offset 0,-0.5
#set ytics font "Arial,20" offset 0,0

#set tmargin 0
#set bmargin 0
#set lmargin 5
#set rmargin 0
unset xtics
unset ytics

set multiplot layout 2, 2

set key autotitle column nobox samplen 1 noenhanced
unset title

set pm3d map
set palette color negative

# [1,1]
set tmargin at screen 0.95; set bmargin at screen 0.55
set lmargin at screen 0.08; set rmargin at screen 0.46
set label "elastico" at graph 0.42, graph 1.02 font "Arial Black,25"
unset cblabel
unset colorbox
set ytics nomirror font "Arial,20" offset 0,0
set ylabel "V (Ry)" offset 0,0 font "Arial,25"
splot fileEl
unset label

# [1,2]
set tmargin at screen 0.95; set bmargin at screen 0.55
set lmargin at screen 0.48; set rmargin at screen 0.86
set label "total" at graph 0.44, graph 1.02 font "Arial Black,25"
set cblabel labelCB offset 3.5,0 font "Arial,25"
set colorbox
unset ytics
unset ylabel
splot fileTot
unset label

# [2,1]
set tmargin at screen 0.50; set bmargin at screen 0.10
set lmargin at screen 0.08; set rmargin at screen 0.46
set label "inel simetrico" at graph 0.36, graph 1.02 font "Arial Black,25"
unset cblabel
unset colorbox
set xlabel "E_F (Ry)" font "Arial,25"
set xtics nomirror font "Arial,20" offset 0,0
set ytics nomirror font "Arial,20" offset 0,0
set ylabel "V (Ry)" offset 0,0 font "Arial,25"
splot fileSy
unset label

# [2,2]
set tmargin at screen 0.50; set bmargin at screen 0.10
set lmargin at screen 0.48; set rmargin at screen 0.86
set label "inel assimetrico" at graph 0.35, graph 1.02 font "Arial Black,25"
set cblabel labelCB offset 0,0 font "Arial,25"
set colorbox
unset ytics
unset ylabel
set xlabel "E_F (Ry)" font "Arial,25"
splot fileAsy

unset multiplot

