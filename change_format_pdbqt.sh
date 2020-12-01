#!/bin/bash

for file in ./*
do
    if [ ${file#*.} == 'mol2' ]
    then  
        NAME=${file%%.*}
        obabel -imol2 $file -opdbqt $NAME.pdbqt
        rm $file
    fi
done  