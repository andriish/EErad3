#!/bin/bash
#for SL7   

CARDS=( E00.y1d8.iL0.LO.input   
        Etx.y1d8.iN1.NLOicol1.input  Etx.y1d8.iN2.NLOicol2.input Etx.y1d8.iN3.NLOicol3.input  
        Etx.y1d6.iZ1.NNLOicol1.input  
        Etx.y1d6.iZ2.NNLOicol2.input  
        Etx.y1d6.iZ3.NNLOicol3.input  
        Etx.y1d6.iZ4.NNLOicol4.input  
        Etx.y1d6.iZ5.NNLOicol5.input 
        Etx.y1d6.iZ6.NNLOicol6.input  
        )
mkdir temp
Q=90
for a in $(ls calc/E*.y1d6.iZ1.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ1.$Q.txt
for a in $(ls calc/E*.y1d6.iZ2.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ2.$Q.txt
for a in $(ls calc/E*.y1d6.iZ3.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ3.$Q.txt
for a in $(ls calc/E*.y1d6.iZ4.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ4.$Q.txt
for a in $(ls calc/E*.y1d6.iZ5.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ5.$Q.txt
for a in $(ls calc/E*.y1d6.iZ6.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iZ6.$Q.txt


for a in $(ls calc/E*.y1d8.iN1.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iN1.$Q.txt
for a in $(ls calc/E*.y1d8.iN2.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iN2.$Q.txt
for a in $(ls calc/E*.y1d8.iN3.ELa); do  echo -n $(echo $a | cut -c7-10); sed "${Q}q;d" $a; done  > temp/iN3.$Q.txt
