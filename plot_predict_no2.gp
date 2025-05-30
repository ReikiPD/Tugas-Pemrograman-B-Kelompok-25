set datafile separator ","
set title 'Prediksi Konsentrasi NO_{2} terhadap Waktu'
set xlabel 'Waktu (jam)'
set ylabel 'NO2 (ppb)'
set grid
set xrange [8:16]
plot \
  'predict_no2_data.csv' using 1:2 with lines lw 2 title 'Lagrange', \
  ''                        using 1:3 with lines lw 2 title 'Spline Kubik', \
  ''                        using 1:4 with lines lw 2 title 'Regresi Kuadrat'
pause -1
