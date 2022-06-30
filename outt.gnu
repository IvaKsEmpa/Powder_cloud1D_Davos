b = 1000
reset
set terminal postscript enhanced eps color font 'Helvetica, 12.5'
#set term png medium size 640,480 
set encoding utf8
set style line 5 lt rgb "red" lw 1 pt 1
#set term png  font arial 14 size 800,600
set style data linespoints
#set terminal postscript enhanced 
set border linewidth 1.5
set grid
#set parametric
set nokey
set key left # bottom 
#set xrange[1500.:2000.]
#set yrange[0.0:0.1]  
v= sprintf("%g",a)
set output './images/powder628_hHath'.v.'.eps'   
#set multiplot layout 4,1 columns
set xlabel 'x (m)'
set ylabel "h (m)"
plot './resu/powder628.txt' u 1:6  every:::a::a  lt rgb 'red' lw 2 with lines  title 'Numerical height hat h', './resu/powder628.txt' u 1:2  every:::a::a  lt rgb 'blue' lw 2 with lines  title 'Numerical height h'
set ylabel "U (m/s)"
#plot './resu/powder509TimeEvol.txt' u 1:3 every:::a::a  lt -1  lw 2 title 'Numerical velocity' with lines
set ylabel "kinetic ener"
#plot './resu/powder509TimeEvol.txt' u 1:(($3*$3/2)) every:::a::a  lt -1 lw 2 title 'Numerical kinetic energy' with lines 
set ylabel "turb ener"
#plot './resu/powder509TimeEvol.txt' u 1:((2+$4)*$2*$2/2) every:::a::a  lt -1  lw 2 title 'Numerical turbulent energy'   with lines 
#set ylabel "pulse"
#plot './resu/powder509TimeEvol.txt' u 1:5 every:::a::a  lt -1 lw 2 title 'Pulse function'   with lines 
#set ylabel "Q"
#plot './resu/qdt.out' u 1:2 every:::a::a  lt -1  lw 2 title 'Q'   with lines 
#unset multiplot
pause -1
a=a+1
if(a<b) reread


 # 'LineProfiles.TXT' u ($13+1200):14, 'LineProfiles.TXT' u ($23 +1300):24,  'LineProfiles.TXT' u ($9 +1070):10, 'LineProfiles.TXT' u ($11+1150):12,
 #plot 'LineProfiles.TXT' u ($23 +1300):24,  'LineProfiles.TXT' u ($9 +1070):10, 'LineProfiles.TXT' u ($11+1150):12,'LineProfiles.TXT' u ($13+1200):14,  #PulseSin.out' u 1:6 every:::a::a lt rgb 'red'  lw 2  

#set ylabel "kinetic energy"
#plot './resu/test1DCr1.out' u 1:(($2*$3*$3/2)) every:::a::a  lt rgb 'green'  title 'Numerical results' 
#set ylabel "turbulent energy"
#plot './resu/test1DCr1.out' u 1:($4*$2*$2*$2/2) every:::a::a  lt rgb 'green'  title 'Numerical results'  




#plot './resu/PulseSin.out' u 1:5 every:::a::a with lines title 'core'
set ylabel "h (m)"
#plot 'LineProfiles.TXT' u ($19+1280):20 with lines lt rgb 'blue' title 'ground-based photogrammetry measurements', './resu/PulseSin.out' u 1:6 every:::a::a with lines lw 2 lt rgb 'black' title 'numerical', 'LineProfiles.TXT' u ($23 +1300):24 with lines lt rgb 'red'
set ylabel "U (m/s)"
# plot './resu/PulseSin.out' u 1:3 every:::a::a  lt rgb 'red'  lw 2 title 'numerical'with lines
set ylabel "Phi"
# plot './resu/PulseSin.out' u 1:4 every:::a::a  lt rgb 'blue'  title 'numerical' with lines 
unset multiplot
pause -1
a=a+1
if(a<b) reread

#set ylabel "h (m)"
#plot 'LineProfiles.TXT' u ($13+1270):14 pt 2 title 'Measurements',  './resu/TestEnergy.out' u 1:6 every:::a::a with lines lw 2 lt -1 title 'numerical height' #,  'LineProfiles.TXT' u ($9 +1120):10, 'LineProfiles.TXT' u ($11+1180):12 title '45 s'
#set ylabel "U (m/s)"
#plot './resu/TestEnergy.out' u 1:3 every:::a::a  lt -1  lw 2 title 'Numerical velocity' with lines
set ylabel "Phi (1/s^{2})"
#plot './resu/TestEnergy.out' u 1:4 every:::a::a  lt -1  title 'Numerical results' with lines 
set ylabel "kin energy (m^{2}/s^{2})"
#plot './resu/TestEnergy.out' u 1:(($3*$3/2)) every:::a::a  lt -1  title 'Numerical kinetic energy'  with lines 
set ylabel "turb energy (m^{2}/s^{2})"
#plot './resu/TestEnergy.out' u 1:($4*$2*$2/2) every:::a::a  lt -1  title 'Numerical turbulent energy'   with lines 

set ylabel "h (m)"
#plot '509_alongtrack_profiles.txt' u ($1+1230):2 pt 2 title 'Measurements', './resu/Test509.out' u 1:6 every:::a::a with lines lw 2 lt -1 title 'numerical height' #,  'LineProfiles.TXT' u ($9 +1120):10, 'LineProfiles.TXT' u ($11+1180):12 title '45 s'
set ylabel "U (m/s)"
#plot './resu/Test509.out' u 1:3 every:::a::a  lt -1  lw 2 title 'Numerical velocity' with lines
set ylabel "Phi (1/s^{2})"
#plot './resu/Test509.out' u 1:4 every:::a::a  lt -1  title 'Numerical results' with lines 
set ylabel "kin energy (m^{2}/s^{2})"
#plot './resu/Test509.out' u 1:(($3*$3/2)) every:::a::a  lt -1  title 'Numerical kinetic energy'  with lines 
set ylabel "turb energy (m^{2}/s^{2})"
#plot './resu/Test509.out' u 1:($4*$2*$2/2) every:::a::a  lt -1  title 'Numerical turbulent energy'   with lines 