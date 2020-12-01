#! /bin/bash

inputdir=$1
confile=$2
outdir=$3

echo $inputdir
#read -p "input work dir and config file" inputdir confile 
#read -p "input output dir" outdir

for f in $inputdir/*.pdbqt; do
    b=`basename $f .pdbqt`
    # echo Processing ligand $b
    if [ -d $outdir/$b ] && [ `ls $outdir/$b | wc -w` = 2 ];then 
       continue
    fi
    echo 'vina_screen.sh mkdir -p $$outdir/$b' $outdir/$b
    mkdir -p $outdir/$b
    vina --config $confile --ligand $f --out $outdir/${b}/out.pdbqt --log $outdir/${b}/log.txt
done
