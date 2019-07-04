#!/bin/bash
set -x
ZMCSP=EERAD3000.tgz

tar xvzf $ZMCSP
ls -lah
pwd
uname -a
ldconfig -v
find 

CARD=$(cat dir.txt | head -n 1| tail -n 1)
SEED=$(cat dir.txt | head -n 2| tail -n 1)
TOP=$(pwd)
cd $TOP


cd GP

mkdir -p calc
mkdir -p log
cd calc
../eerad3 -i ../cards/$CARD -n $SEED ../log/$CARD'_'$SEED'.log'

cd ..

tar cfz $CARD'_'$SEED'.tgz' log calc 

mv $CARD'_'$SEED'.tgz' $TOP
