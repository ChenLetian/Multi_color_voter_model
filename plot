set style line 1 lc rgb '#D53E4F' lt 1 lw 10 # red
set style line 2 lc rgb '#F46D43' lt 1 lw 10 # orange
set style line 3 lc rgb '#FDAE61' lt 1 lw 10 # pale orange
set style line 4 lc rgb '#FEE08B' lt 1 lw 10 # pale yellow-orange
set style line 5 lc rgb '#E6F598' lt 1 lw 10 # pale yellow-green
set style line 6 lc rgb '#ABDDA4' lt 1 lw 10 # pale green
set style line 7 lc rgb '#66C2A5' lt 1 lw 10 # green
set style line 8 lc rgb '#3288BD' lt 1 lw 10 # blue
set term postscript enhanced color
set xlabel 'Time'
set ylabel 'Probabilities'
set output '3_species.eps'
plot 'Collision_100000_3_2_5000_3.txt' u 1:($2/$4) w l ls 2 title 'P(annhilation|collision)',\
'Collision_100000_3_2_5000_3.txt' u 1:($3/$4) w l ls 8 title 'P(coalescence|collision)',\
'Collision_100000_3_2_5000_3.txt' u 1:(($2-$3)/$4) w l ls 7 title 'Delta'
set output '4_species.eps'
plot 'Collision_100000_4_2_5000_3.txt' u 1:($2/$4) w l ls 2 title 'P(annhilation|collision)',\
'Collision_100000_4_2_5000_3.txt' u 1:($3/$4) w l ls 8 title 'P(coalescence|collision)',\
'Collision_100000_4_2_5000_3.txt' u 1:(($2-$3)/$4) w l ls 7 title 'Delta'
set output '5_species.eps'
plot 'Collision_100000_5_2_5000_3.txt' u 1:($2/$4) w l ls 2 title 'P(annhilation|collision)',\
'Collision_100000_5_2_5000_3.txt' u 1:($3/$4) w l ls 8 title 'P(coalescence|collision)',\
'Collision_100000_5_2_5000_3.txt' u 1:(($2-$3)/$4) w l ls 7 title 'Delta'

set log
set output '3_species_triplets.eps'
plot 'Collision_100000_3_2_5000_3.txt' u 1:($5/$8) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_3_2_5000_3.txt' u 1:($6/$8) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'
set output '4_species_triplets.eps'
plot 'Collision_100000_4_2_5000_3.txt' u 1:($5/$8) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_4_2_5000_3.txt' u 1:($6/$8) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'
set output '5_species_triplets.eps'
plot 'Collision_100000_5_2_5000_3.txt' u 1:($5/$8) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_5_2_5000_3.txt' u 1:($6/$8) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'

set output '3_species_triplets_relative.eps'
plot 'Collision_100000_3_2_5000_3.txt' u 1:($5/$7) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_3_2_5000_3.txt' u 1:($6/$7) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'
set output '4_species_triplets_relative.eps'
plot 'Collision_100000_4_2_5000_3.txt' u 1:($5/$7) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_4_2_5000_3.txt' u 1:($6/$7) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'
set output '5_species_triplets_relative.eps'
plot 'Collision_100000_5_2_5000_3.txt' u 1:($5/$7) w l ls 2 title 'P(iji)/(P(iji)+P(ijk))',\
'Collision_100000_5_2_5000_3.txt' u 1:($6/$7) w l ls 7 title 'P(ijk)/(P(iji)+P(ijk))'


//New gnuplot!!!
set style line 1  lc rgb '#0b1d6a' lt 1 lw 5 
set style line 2  lc rgb '#0e2588' lt 1 lw 5
set style line 3  lc rgb '#102b9f' lt 1 lw 5
set style line 4  lc rgb '#1437cd' lt 1 lw 5
set style line 5  lc rgb '#163ce1' lt 1 lw 5
set style line 6  lc rgb '#1641f8' lt 1 lw 5
set style line 7  lc rgb '#1652f8' lt 1 lw 5
set style line 8 lc rgb '#1662f8' lt 1 lw 5
set style line 9 lc rgb '#1672f8' lt 1 lw 5
set style line 10 lc rgb '#1682f8' lt 1 lw 5
set style line 11  lc rgb '#1692f8' lt 1 lw 5
set style line 12  lc rgb '#169ef8' lt 1 lw 5
set style line 13  lc rgb '#16b1f8' lt 1 lw 5
set style line 14  lc rgb '#16c3f8' lt 1 lw 5
set style line 15  lc rgb '#16cff8' lt 1 lw 5
set style line 16  lc rgb '#16dbf8' lt 1 lw 5
set style line 17  lc rgb '#16ecf8' lt 1 lw 5
set style line 18  lc rgb '#16f8f8' lt 1 lw 5
set style line 19  lc rgb '#16f8df' lt 1 lw 5
set style line 20 lc rgb '#16f8ca' lt 1 lw 5
set style line 21 lc rgb '#16f8b6' lt 1 lw 5
set style line 22 lc rgb '#16f8a1' lt 1 lw 5
set style line 23  lc rgb '#16f886' lt 1 lw 5
set style line 24  lc rgb '#16f86f' lt 1 lw 5
set style line 25  lc rgb '#16f855' lt 1 lw 5
set style line 26  lc rgb '#16f839' lt 1 lw 5
set style line 27  lc rgb '#16f823' lt 1 lw 5
set style line 28  lc rgb '#16f811' lt 1 lw 5
set log
set term postscript enhanced color
unset key
set output '3_species_detailled.eps'
plot for [col=2:13] 'Collision_detailed_100000_3_2_5000_4.txt' using 1:(column(col)/($14)) with lines lw 5
set output '4_species_detailled.eps'
plot for [col=2:13] 'Collision_detailed_100000_4_2_5000_4.txt' using 1:(column(col)/($38)) with lines lw 5
set output '5_species_detailled.eps'
plot for [col=2:13] 'Collision_detailed_100000_5_2_5000_4.txt' using 1:(column(col)/($82)) with lines lw 5











