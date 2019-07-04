#!/bin/bash
b=$(basename $1)
../eerad3 -i $1 -n $2 &> ../log/$b'_'$2'.log' 