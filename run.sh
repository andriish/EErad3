#!/bin/bash
set -x
mkdir -p calc
mkdir -p log
mkdir -p logcombine
cd calc

for a in $(ls -1 ../cards/E0*input| grep ZZZZ ); do
:
b=$(basename $a)
for i in `seq 200 205`; do
../eerad3 -i $a -n $i &> ../log/$b'_'$i'.log' &
done
done


for a in $(ls -1 ../cards/Et*input | grep '\.NLO'| grep ZZZZ); do
:
b=$(basename $a)
for i in `seq 200 205`; do
../eerad3 -i $a -n $i &> ../log/$b'_'$i'.log' &
done
done

for a in $(ls -1 ../cards/Et*input | grep '\.NNLO' | grep ZZZ); do
:
b=$(basename $a)
for i in `seq 200 205`; do
../eerad3 -i $a -n $i &> ../log/$b'_'$i'.log' &
done
done


wait 
for a in $(ls -1 ../cards/E*combine); do
b=$(basename $a)
../eerad3_combine -i $a  &> ../logcombine/$b'.log'
done

../eerad3_dist -i ../cards/eerad3_dist.input
#cd doc
#sh bin/convert_predictions.sh