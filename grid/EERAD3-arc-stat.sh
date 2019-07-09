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
        

for CARD in "${CARDS[@]}"
do

echo $CARD $(ls output/$CARD* | wc -l)
done
