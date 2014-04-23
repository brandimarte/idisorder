# --- Input file ---
# Pass 'file' in command line as: gnuplot -e "file='ExVxIsy.tot'" graph3D.gnu > fig.jpg
#file='ExVxIsy.tot'

#E_F = -4.335492603 eV
#Eph =  0.454431528 eV

# --- Title ---
#set key autotitle column nobox samplen 1 noenhanced
unset title
#set title "Elastic\n{/*0.6E_F  = -0.0234 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
#set title "Total (elastic + inelastic)\n{/*0.6E_F  = -0.0234 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
#set title "Inelastic (symmetric part)\n{/*0.6E_F  = -0.0234 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
#set title "Inelastic (asymmetric part)\n{/*0.6E_F  = -0.0234 eV}\n{/*0.6E_{ph} =  0.4544 eV}"

# --- Size and output characteristics ---
set terminal jpeg transparent enhanced font "Arial,45" fontscale 1.0 size 4400, 3600

# --- 3D position ---
set view 70, 20, 1, 1
#set view 67, 22, 1, 1

# --- Axis ranges ---
set xrange [ -1.0000 : 1.0000 ] noreverse nowriteback
#set xrange [ 0.000 : 1.80000 ] noreverse nowriteback
set yrange [ -0.50 : 0.50 ] noreverse nowriteback
set zrange [ -2.0 : 3.0 ] noreverse nowriteback
set cbrange [ -4.0 : 4.0 ]

# --- Tics ---
set xtics font "Arial,35" offset 0,-0.5
set ytics font "Arial,35" offset -0.1,-0.4
set ztics font "Arial,35" offset 0.2,0
#set ztics -1.00000,0.25,1.00000 norangelimit
#unset xtics
#unset ytics
set xtics out nomirror
set ytics out nomirror
set ztics out
set format x "%2.1f"
set format y "%2.1f"

# --- Grid lines ---
set grid ztics lw 2
set ticslevel 0
#set dgrid3d 30,30
#set isosamples 10,10

# --- Borders ---
#set border 1+2+4+8+16+32+64+256+512
set border 1+2+4+8+16+32+64

# --- Margin sizes ---
set tmargin 0
set bmargin 0
set lmargin 5
set rmargin 0

# --- Multiple plots ---
#set multiplot layout 2, 2

# --- Plot style ---
set pm3d at b
#set palette defined (-0.004 "#696969", 0 "#F5F5F5", 0.004 "#DC143C")
set palette color negative

# --- Axis labels ---
#set tmargin at screen 0.95; set bmargin at screen 0.55
#set lmargin at screen 0.08; set rmargin at screen 0.46
unset clabel # remove nome da legenda de cores
unset colorbox # remove legenda de cores
unset key # remove legenda
set xlabel "E - E_F (eV)" 
set xlabel  offset character -2, -2, 0 font "Arial,45" textcolor lt -1 norotate
set ylabel "V (V)" 
set ylabel  offset character 1, -1, 0 font "Arial,45" textcolor lt -1 rotate by -270
set zlabel "dI/dV (G__0)" 
set zlabel  offset character -2, 0, 0 font "Arial,45" textcolor lt -1 rotate by 90

# --- Plot ---
#splot file with lines lt 2 lw 0
#splot file with linespoints ls 2

# --- Multiplot ---
set multiplot layout 2, 2
#fileEl='ExVxdIel.tot'
#fileTot='ExVxdItot.tot'
#fileSy='ExVxdIsy.tot'
#fileAsy='ExVxdIasy.tot'
#fileEl='ExVxdIel.red1'
#fileTot='ExVxdItot.red1'
#fileSy='ExVxdIsy.red1'
#fileAsy='ExVxdIasy.red1'
#fileEl='ExVxdIel.red2'
#fileTot='ExVxdItot.red2'
#fileSy='ExVxdIsy.red2'
#fileAsy='ExVxdIasy.red2'
fileEl='ExVxdIel.red3'
fileTot='ExVxdItot.red3'
fileSy='ExVxdIsy.red3'
fileAsy='ExVxdIasy.red3'
# [1,1]
#set title "Elastic\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleEl
splot fileEl with lines lt 2 lw 0
# [1,2]
#set title "Total (elastic + inelastic)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleTot
splot fileTot with lines lt 2 lw 0
# [2,1]
set zrange [ -0.001 : 0.001 ] noreverse nowriteback
set cbrange [ -3.0e-04 : 3.0e-04 ]
set format z "%.0e"
#set title "Inelastic (symmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleIsy
splot fileSy with lines lt 2 lw 0
# [2,2]
#set title "Inelastic (asymmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleIasy
splot fileAsy with lines lt 2 lw 0
