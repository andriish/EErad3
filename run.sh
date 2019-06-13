#!/bin/bash
for a in $(ls -1 cards/E0*); do
./eerad3 -i $a -n 0 &
done


#!/bin/bash
for a in $(ls -1 cards/Et*); do
./eerad3 -i $a -n 1 &
./eerad3 -i $a -n 2 &
./eerad3 -i $a -n 3 &
done


wait 


./eerad3_combine



./eerad3_dist