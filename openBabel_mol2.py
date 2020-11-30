#!/usr/local/bin/python3
import os
import subprocess
import sys 
import shutil
from config import Config

filePath=Config.DIR_PATH
dir_list=os.listdir(filePath)

def call_perl(matchpath,file):
    print('start running perl script...')
    ret=subprocess.call(['perl',matchpath+'/splitmol2files.pl',file])
    if ret==0:
        print('success')
    else:
        print('failed')

def change_format(matchPath):
    print('start running shell script')
    os.system(matchPath+'/change_format_pdbqt.sh')
    # ret=subprocess.call([matchPath+'/change_format_pdbqt.sh'])
    # if ret==0:
    #     print('success')
    # else:
    #     print('failed')

def one_2_infinite(chdir,matchpath,file_list):
    for file in file_list:
        if file == 'actives_final.mol2':
            if not os.path.exists(chdir+'/ligand'):
                os.mkdir(chdir+'/ligand')
            os.chdir(chdir+'/ligand')
            shutil.move(chdir+'/'+file,chdir+'/ligand/'+file)
        else:
            if not os.path.exists(chdir+'/decoy'):
                os.mkdir(chdir+'/decoy')
            os.chdir(chdir+'/decoy')
            shutil.move(chdir+'/'+file,chdir+'/decoy/'+file)
        call_perl(matchpath,file)
        change_format(matchpath)



for i in dir_list:
    if not os.path.isdir(filePath+'/'+i):
        continue
    one_2_infinite(filePath+'/'+i+'/dude','/root/lcq/transform_mol2_pdbqt',['actives_final.mol2','decoys_final.mol2'])
