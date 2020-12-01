#! /bin/bash

inputdir=$1
# input example:/home/leicq/flexbenchmark/data/gpcr
wkdir=$(pwd)
function dockingFuc(){
    $wkdir/vina_screen_local.sh $1 $2 $3 &
}


# change whole protein dock input in gpcr dir into pdbqt format 
#echo 'dock_loop.sh: stage1 change format'
#python3 ./openBabel_mol2.py

echo 'dock_loop.sh: stage2 vina virtual docking' 
i=0
for protein in $inputdir/*;do
    proteinclass[i]=$(basename $protein)
    i='expr $i + 1'
done
#proteinclass=`for f in $inputdir/*; do bn=$(basename $f);done | uniq`

echo protein_class is $proteinclass
for pr in $proteinclass; do
    echo do $pr Vs 
    dockingFuc $inputdir/$pr/dock/input/ligand $inputdir/$pr/dock/input/rigid_conf.txt $inputdir/$pr/dock/output/rigid/ligand
   # dockingFuc $inputdir/$pr/dock/input/decoy $inputdir/$pr/dock/input/rigid_conf.txt $inputdir/$pr/dock/output/rigid/decoy
    dockingFuc $inputdir/$pr/dock/input/ligand $inputdir/$pr/dock/input/flex_conf.txt $inputdir/$pr/dock/output/flex/ligand 
    wait
done

for pr in $proteinclass; do
    dockingFuc $inputdir/$pr/dock/input/decoy $inputdir/$pr/dock/input/rigid_conf.txt $inputdir/$pr/dock/output/rigid/decoy
    dockingFuc $inputdir/$pr/dock/input/decoy $inputdir/$pr/dock/input/flex_conf.txt $inputdir/$pr/dock/output/flex/decoy 
    wait
    echo 'dock_loop.sh: stage3 paint roc curve'
    python3 $wkdir/generatedockres.py $inputdir/$pr
done



