#!/bin/bash
set -x
mkdir -p calc
mkdir -p log
mkdir -p logcombine
cd calc
#for a in $(ls -1 ../cards/E0*input | grep QQQQ); do
for a in $(ls -1 ../cards/E0*input ); do
:
b=$(basename $a)
../eerad3 -i $a -n 0 &> ../log/$b'_0.log' & 
#../eerad3 -i $a -n 1 &> ../log/$b'_1.log' & 
done

#wait
#exit

#!/bin/bash

for a in $(ls -1 ../cards/Et*input | grep '\.NLO'); do
#for a in $(ls -1 ../cards/Et*input | grep '\.NLO' | grep QQQQ ); do
:
b=$(basename $a)
../eerad3 -i $a -n 1 &> ../log/$b'_1.log' &
../eerad3 -i $a -n 2 &> ../log/$b'_2.log' &
../eerad3 -i $a -n 3 &> ../log/$b'_3.log' &
../eerad3 -i $a -n 4 &> ../log/$b'_4.log' &
../eerad3 -i $a -n 5 &> ../log/$b'_5.log' &
../eerad3 -i $a -n 6 &> ../log/$b'_6.log' &
../eerad3 -i $a -n 7 &> ../log/$b'_7.log' &
../eerad3 -i $a -n 8 &> ../log/$b'_8.log' &
../eerad3 -i $a -n 9 &> ../log/$b'_9.log' &
../eerad3 -i $a -n 10 &> ../log/$b'_10.log' &
../eerad3 -i $a -n 11 &> ../log/$b'_11.log' &
../eerad3 -i $a -n 12 &> ../log/$b'_12.log' &
../eerad3 -i $a -n 13 &> ../log/$b'_13.log' &
../eerad3 -i $a -n 14 &> ../log/$b'_14.log' &
../eerad3 -i $a -n 15 &> ../log/$b'_15.log' &
done

#for a in $(ls -1 ../cards/Et*input | grep '\.NNLO'| grep QQQQ); do
for a in $(ls -1 ../cards/Et*input | grep '\.NNLO'); do
:
b=$(basename $a)
../eerad3 -i $a -n 1 &> ../log/$b'_1.log' &
../eerad3 -i $a -n 2 &> ../log/$b'_2.log' &
../eerad3 -i $a -n 3 &> ../log/$b'_3.log' &
../eerad3 -i $a -n 4 &> ../log/$b'_4.log' &
../eerad3 -i $a -n 5 &> ../log/$b'_5.log' &
../eerad3 -i $a -n 6 &> ../log/$b'_6.log' &
../eerad3 -i $a -n 7 &> ../log/$b'_7.log' &
../eerad3 -i $a -n 8 &> ../log/$b'_8.log' &
../eerad3 -i $a -n 9 &> ../log/$b'_9.log' &
../eerad3 -i $a -n 10 &> ../log/$b'_10.log' &
../eerad3 -i $a -n 11 &> ../log/$b'_11.log' &
../eerad3 -i $a -n 12 &> ../log/$b'_12.log' &
../eerad3 -i $a -n 13 &> ../log/$b'_13.log' &
../eerad3 -i $a -n 14 &> ../log/$b'_14.log' &
../eerad3 -i $a -n 15 &> ../log/$b'_15.log' &
done


wait 
for a in $(ls -1 ../cards/Et*combine); do
b=$(basename $a)
../eerad3_combine -i $a  &> ../logcombine/$b'.log'
done

../eerad3_dist -i ../cards/eerad3_dist.input