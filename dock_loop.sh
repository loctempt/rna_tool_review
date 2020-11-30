#!bin/bash

inputdir=$1
# input example:/home/leicq/benchmark/data/gpcr

dockingFuc(inputdir,confile,outdir){
    fork ./vina_screen_local.sh $indir $confile $outdir
}

obableFuc(){
    
}

proteinclass=`for f in $inputdir/*; do bn=$(basename $f); echo ${bn}; done | uniq`

for pr in $proteinclass; do
    


