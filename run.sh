#!/bin/bash
for a in $(ls -1 cards/E0*input); do
./eerad3 -i $a -n 0 &
done


#!/bin/bash
for a in $(ls -1 cards/Et*input); do
./eerad3 -i $a -n 1 &
./eerad3 -i $a -n 2 &
./eerad3 -i $a -n 3 &
done


wait 
for a in $(ls -1 cards/Et*combine); do
./eerad3_combine -i $a 
done

./eerad3_dist