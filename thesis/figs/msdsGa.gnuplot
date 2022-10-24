set datafile separator ','

set xlabel "Generation"

plot "msdsGa.csv" using 1:2 title "Sum of fittest population loss" with lines