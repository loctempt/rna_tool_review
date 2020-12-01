#!/bin/bash
activedir=$1
decoydir=$2

function change_format(inputdir){
   cd $inputdir
   for file in $inputdir/*;do
    if [ ${file#*.} == 'mol2' ]
    then  
        NAME=${file%%.*}
        obabel -imol2 $file -opdbqt -O $NAME.pdbqt
        rm $file
    fi
    done  
}

change_format($activedir)
change_format($decoydir)
