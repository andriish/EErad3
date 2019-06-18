#!/bin/bash
set -x
mkdir -p calc
mkdir -p log
mkdir -p logcombine
cd calc
for a in $(ls -1 ../cards/E0*input); do
:
b=$(basename $a)
../eerad3 -i $a -n 0 &> ../log/$b'_0.log' & 
#../eerad3 -i $a -n 1 &> ../log/$b'_1.log' & 
done

#wait
#exit

#!/bin/bash
for a in $(ls -1 ../cards/Et*input | grep '\.NLO'); do
#for a in $(ls -1 ../cards/Et*input | grep '\.NNLO'); do
#for a in $(ls -1 ../cards/Et*input ); do
:
b=$(basename $a)
../eerad3 -i $a -n 1 &> ../log/$b'_1.log' &
../eerad3 -i $a -n 2 &> ../log/$b'_2.log' &
../eerad3 -i $a -n 3 &> ../log/$b'_3.log' &
../eerad3 -i $a -n 4 &> ../log/$b'_4.log' &
../eerad3 -i $a -n 5 &> ../log/$b'_5.log' &
done


wait 
for a in $(ls -1 ../cards/Et*combine); do
../eerad3_combine -i $a  &> ../logcombine/$a'.log'
done

../eerad3_dist -i ../cards/eerad3_dist.input