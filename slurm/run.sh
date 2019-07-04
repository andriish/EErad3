#!/bin/bash
set -x
mkdir -p calc
mkdir -p log
mkdir -p logcombine
cd calc
cp ../runner.sh ./
set -x
#echo '
#!#/bin/bash
#b=$(basename $1)
#../eerad3 -i $1 -n $2 &> ../log/$b'_'$2'.log' 
#'
# > runner.sh
chmod +x runner.sh

for a in $(ls -1 ../cards/E0*input| grep QQ ); do
:
b=$(basename $a)
for i in `seq 1 200`; do
:
sbatch -p standard --workdir=$(pwd)  runner.sh  $a $i
done
done
#exit
for a in $(ls -1 ../cards/Et*input | grep '\.NLO'); do
:
b=$(basename $a)
for i in `seq 1 100`; do
sbatch -p standard --workdir=$(pwd)  runner.sh  $a $i
done
done

for a in $(ls -1 ../cards/Et*input | grep '\.NNLO'| grep QQ); do
:
for i in `seq 1 50`; do
:
sbatch -p standard --workdir=$(pwd)  runner.sh  $a $i
done
done
