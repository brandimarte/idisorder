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
#set view 90, 90, 1, 1 # eixo de bias
#set view 90, 0, 1, 1 # eixo de energia
#set view 67, 22, 1, 1

# --- Axis ranges ---
set xrange [ -2.0000 : 2.0000 ] noreverse nowriteback
#set xrange [ -1.1 : -1.0 ] noreverse nowriteback
set yrange [ -0.50 : 0.50 ] noreverse nowriteback
set zrange [ -4.0 : 4.0 ] noreverse nowriteback
set cbrange [ -0.2 : 0.2 ]

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
set zlabel "(d^2I/dV^2)/(dI/dV) (1/V)" 
set zlabel  offset character -2, 0, 0 font "Arial,45" textcolor lt -1 rotate by 90

# --- Plot ---
#splot file with lines lt 2 lw 0
#splot file with linespoints ls 2

# --- Multiplot ---
set multiplot layout 2, 2
#fileEl='IETSel.tot'
#fileTot='IETStot.tot'
#fileSy='IETSsy.tot'
#fileAsy='IETSasy.tot'
#fileEl='IETSel.red1'
#fileTot='IETStot.red1'
#fileSy='IETSsy.red1'
#fileAsy='IETSasy.red1'
#fileEl='IETSel.red2'
#fileTot='IETStot.red2'
#fileSy='IETSsy.red2'
#fileAsy='IETSasy.red2'
fileEl='IETSel.red3'
fileTot='IETStot.red3'
fileSy='IETSsy.red3'
fileAsy='IETSasy.red3'
# [1,1]
#set cbrange [-1.0:1.0]
#set title "Elastic\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleEl
splot fileEl with lines lt 2 lw 0
# [1,2]
#set title "Total (elastic + inelastic)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleTot
splot fileTot with lines lt 2 lw 0
# [2,1]
set zrange [ -600.0 : 800.0 ] noreverse nowriteback
set cbrange [ -500.0 : 500.0]
#set title "Inelastic (symmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleIsy
splot fileSy with lines lt 2 lw 0
# [2,2]
set zrange [ -6000.0 : 6000.0 ] noreverse nowriteback
set cbrange [ -4000.0 : 4000.0]
#set title "Inelastic (asymmetric part)\n{/*0.6E_F  = -4.3354 eV}\n{/*0.6E_{ph} =  0.4544 eV}"
set title titleIasy
splot fileAsy with lines lt 2 lw 0
