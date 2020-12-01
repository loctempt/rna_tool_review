#!/usr/local/bin/python3
import os
import subprocess
import sys 
import shutil
import gzip
from config import Config

filePath=Config.DIR_PATH
dir_list=os.listdir(filePath)
 #TODO 添加判断小分子pdbqt中是否有不可以执行的小分子，并剔除；

def ugzip(filename):
    file=filename.replace('.gz','')
    gunzip=gzip.GzipFile(filename)
    open(file, "wb+").write(gunzip.read())
    gunzip.close()

def call_pdbqyt_legalization(path):
    filelist=os.listdir(path)
    for file in filelist:
        del_flag=False
        with open(path+'/'+file,'r') as cur_pdbqt:
            for line in cur_pdbqt.readlines():
                # line=line.split()
                if line[0:4] == 'ATOM' and line[12:14].strip() not in ['H', 'C', 'N', 'O', 'F', 'MG', 'P', 'S', 'CL', 'CA', 'MN', 'FE', 'ZN', 'BR' ,'I' ]:
                    del_flag=True
                    break
        if del_flag == True:
            os.remove(path+'/'+file)

def call_perl(file):
    print('obabel_mol2:start running perl script splitmol2files.pl')
    ret=subprocess.call(['perl',Config.SCRIPT_PATH+'/splitmol2files.pl',file])
    if ret==0:
        print('obabel_mol2 splitmol2files.pl:success')
    else:
        print('splitmol2files.pl: failed')

def change_format(lig_path,dec_path):
    print('obabel_mol2:start running shell script change_format_pdbqt.sh')
    os.system(Config.SCRIPT_PATH+'/change_format_pdbqt.sh '+ lig_path + ' ' + dec_path)
    

def one_2_infinite(chdir,dudedir,file_list):
    ligand_path = chdir+'/ligand'
    decoy_path = chdir+'/decoy'
    # if not os.path.exists(ligand_path): 
    #     os.mkdir(ligand_path)
    # if not os.path.exists(decoy_path) :
    #     os.mkdir(decoy_path)
        
    # for file in file_list:    
    #     if file == 'actives_final.mol2':
    #         os.chdir(ligand_path)
    #         shutil.copy(dudedir+'/'+file,ligand_path+'/'+file)
    #         call_perl(ligand_path+'/'+file)
    #     else:
    #         os.chdir(decoy_path)
    #         shutil.copy(dudedir+'/'+file,decoy_path+'/'+file)
    #         call_perl(decoy_path+'/'+file)
        
    # change_format(ligand_path,decoy_path)
    call_pdbqyt_legalization(ligand_path)
    call_pdbqyt_legalization(decoy_path)


for i in dir_list:
    if not os.path.isdir(filePath+'/'+i):
        continue
    dude_path=filePath+'/'+i+'/dude'
    dock_path=filePath+'/'+i+'/dock/input'
    ugzip(dude_path+'/actives_final.mol2.gz')
    ugzip(dude_path+'/decoys_final.mol2.gz')

    # if not os.path.exists(dude_path):
    #     os.mkdir(dude_path)
    one_2_infinite(dock_path,dude_path,['actives_final.mol2','decoys_final.mol2'])
    os.chdir(Config.SCRIPT_PATH)
