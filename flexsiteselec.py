import os
import math
import random
import numpy as np
from rna_tool_new import * 
import util
from config import Config


pr_path=Config.DIR_PATH
pr_list=os.listdir(pr_path)
pr_dict={}
# 'aa2ar':'3EMLA','adrb1':'2VT4B','adrb2':'3NY8A','braf':'3D4QA',
#     'cdk2':'1H00A','csf1r':'3KRJA','cxcr4':'3ODUB',
#     'drd3':'3PBLA','fak1':'3BZ3A','fgfr1':'3C4FA',
#     'jak2':'3LPBA','kit':'3G0EA','kpcb':'2I0EA',
#     'lck':'2OF2A','mapk2':'3M2WA','mk01':'2OJGA',
#     'mp2k1':'3EQHA','plk1':'2OWBA','rock1':'2ETRA',
#     'src':'3EL8B','tgfr1':'3HMMA','vgfr2':'2P2IB','wee1':'3BIZA',
#     'andr':'2AM9A','esr1':'1SJ0A','esr2':'2FSZA',
#     'gcr':'3BQDA','mcr':'2AA2A','ppara':'2P54A',
#     'pparg':'2GTKA','ppard':'2ZNPA','prgr':'3KBAA',
#     'rxra':'1MV9A','thb':'1Q4XA','akt2':'3D0EA'
    
with open(pr_path+'/class_templateID.txt','r') as class_template_ids:
    for line in class_template_ids.readlines():
        line=line.split()
        pr_dict[line[0]]=line[1]
#key:lig NO (str) -> value: [ligand pdb Pos,small molecule (AAcid) , ctr_point(list[x,y,z])] 
# res=[]

def dist_cutoff_cal(micro_atom_list,macro_atom_list,cutoff):
    tom_cnt=0
    for mi_tom in micro_atom_list:
        if not mi_tom.element_compare('H'):
            continue
        for ma_tom in macro_atom_list:
            if not ma_tom.element_compare('H'):
                continue
            if math.sqrt(mi_tom.atom_dist_2(ma_tom))>cutoff:
                continue
            tom_cnt+=1
            break
    return tom_cnt      

def angle(v1, v2, acute):
# v1 is first vector
# v2 is second vector
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    if (acute == True):
        return angle
    else:
        return 2 * np.pi - angle

def point_divide(point,num):
    return [round(point[0]/num,3),round(point[1]/num,3),round(point[2]/num,3)]

def geo_ctr(lig,return_type):
    sc_X,sc_Y,sc_Z,sc_cnt = 0.0,0.0,0.0,0
    bb_X,bb_Y,bb_Z,bb_cnt = 0.0,0.0,0.0,0
    for atm_rec in lig.get_complete_atom():
        if not atm_rec.element_compare('H'):
           continue
        if atm_rec.get_name().strip() in ['C','N','CA','O']:
            bb_X += atm_rec.get_x()
            bb_Y += atm_rec.get_y()
            bb_Z += atm_rec.get_z()
            bb_cnt+=1
        else:
            sc_X += atm_rec.get_x()
            sc_Y += atm_rec.get_y()
            sc_Z += atm_rec.get_z()
            sc_cnt+=1

    if return_type=='ctr':
        return point_divide([(sc_X+bb_X),(sc_Y+bb_Y),(sc_Z+bb_Z)],(sc_cnt+bb_cnt))
    elif return_type=='bb':
        return point_divide([bb_X,bb_Y,bb_Z],bb_cnt)
    else:
        return point_divide([sc_X,sc_Y,sc_Z],sc_cnt)

def cal_dist(pos_a:list,pos_b:list):
    return round(math.sqrt((pos_a[0] - pos_b[0]) ** 2 +\
                          (pos_a[1] - pos_b[1]) ** 2 +\
                        (pos_a[2] - pos_b[2]) ** 2),2)
                    
def point_bfs(point,molecule_dict):
    if len(molecule_dict) == 0 :
        return None
    cur_class,temp_list=[],[]
    temp_list.append(point)
    
    while len(temp_list)>0:     
        small_molecules_list=list(molecule_dict.keys())
        for key_j in small_molecules_list:
            if molecule_dict[key_j][0]==temp_list[0][1][0]:
                del molecule_dict[point[0]]
                continue
            
            if cal_dist(temp_list[0][1][2],molecule_dict[key_j][2])<=4.5:
                temp_list.append([key_j,molecule_dict[key_j]])
                del molecule_dict[key_j]
        cur_class.append(temp_list.pop(0))
    return cur_class
    
for pr in pr_list:
    cur_big_molecule=None
    # if pr in ['output','test_flex_sc']:
    # if pr not in ['akt2']:
    #     continue
    small_molecules,big_molecules={},{}
    if not os.path.isdir(os.path.join(pr_path,str(pr))):
        continue
    temp=pr_dict[str(pr)]
    # temp=input('plz input '+pr +' template:')
    peptide_list=os.listdir(pr_path+'/'+pr+'/superimpose')
    if not os.path.isdir(pr_path+'/'+pr+'/ligand_test'):
        os.mkdir(pr_path+'/'+pr+'/ligand_test')
    # os.mkdir(pr_path+'\\'+pr+'\\ligand_test\\'+pr)
    # random ligand selection
    lig_NO=0
    for peptide in peptide_list:
        with open(pr_path+'/'+pr+'/superimpose/'+peptide,'r') as input_file:
            # deal big molecule
            cur_peptide=PDB(input_file.readlines())
            if str(peptide[0:5])==temp:
                cur_big_molecule=cur_peptide.get_macro_molecule().get_chain(str(peptide[4]))
            print(peptide[0:5])
            big_molecules[str(peptide[0:5])]=cur_peptide.get_macro_molecule().get_complete_chain()
            # deal small molecules
            cur_small_molecules=cur_peptide.get_micro_molecule()
            if cur_small_molecules == None:
                continue
            # caculate in peptide chain level
            # cur_small_molecules has chain_dict len=1 
            # key:chainID (str) -> value:chain (Chain)
            ligand_count=0
            for small_molecule in cur_small_molecules.chain_dict[str(peptide[4])].get_complete_aa():
                ctr_point=geo_ctr(small_molecule,'ctr')
                small_molecules[str(lig_NO)]=[str(peptide[0:5])+'_'+str(ligand_count),small_molecule,ctr_point]
                ligand_count+=1
                with open(pr_path+'/'+pr+'/ligand_test/'+str(lig_NO)+'.pdb','w') as target:
                    target.writelines(str(small_molecule))
                lig_NO+=1
    
    # point select
    n,m=0,0
    res=[[0 for i in range(lig_NO)] for i in range(lig_NO)]    
    while n < lig_NO:
        while m < lig_NO:
            res[n][m]=cal_dist(small_molecules[str(m)][2],small_molecules[str(n)][2])
            m+=1
        n+=1                
        m=0
    
    with open(pr_path+'/'+pr+'/ligand_test/'+pr+'_detail.txt','w') as detail_f,\
        open(pr_path+'/'+pr+'/ligand_test/'+pr+'_dist.txt','w') as dist_f:
        detail_f.write('Lig_NO\tPDB NAME&lig_num\tctr_x\tctr_y\tctr_z\n')
        for i in range(lig_NO):
            detail_f.write(str(i)+'\t'+small_molecules[str(i)][0]+'\t'+str(small_molecules[str(i)][2][0])+'\t'+str(small_molecules[str(i)][2][1])+'\t'+str(small_molecules[str(i)][2][2])+'\n')
            for j in range(len(res[0])):
              dist_f.write('dist between '+str(i)+' '+str(j)+': '+str(res[i][j])+'\n')    
            dist_f.write('\n')
            
    # cutoff is 10
    # do BFS

    class_list=[]
    while True:
        if small_molecules =={}:
            break
        temp_class=point_bfs([list(small_molecules.keys())[0],small_molecules[list(small_molecules.keys())[0]]],small_molecules)
        if temp_class is not None:        
            class_list.append(temp_class)     
        else:
            break
    
    class_list=sorted(class_list,key=lambda item: len(item),reverse=True)
    
    with open(pr_path+'/'+pr+'/ligand_test/'+pr+'_lig_class.txt','w') as lig_file:
        for i in range(len(class_list)):
            lig_file.write('class'+str(i)+'\n')
            for j in range(len(class_list[i])):
                lig_file.write(class_list[i][j][0]+'\t'+class_list[i][j][1][0]+'\n')

    # cal lig geo 
    i=0
    point_list=[]
    len_flag=False
    print(pr)
    while i<2 :
        lig_x,lig_y,lig_z=0.0,0.0,0.0
        if len_flag is True:
            break
        for j in range(len(class_list[i])) :
            lig_x += class_list[i][j][1][2][0]
            lig_y += class_list[i][j][1][2][1]
            lig_z += class_list[i][j][1][2][2]   
        # solve overfitting problem 
        cur_point=point_divide([lig_x,lig_y,lig_z],len(class_list[i]))
        if len(point_list)==1 and cal_dist(point_list[0],cur_point)<10:
            ref_point=point_divide([point_list[0][0]+cur_point[0],point_list[0][1]+cur_point[1],point_list[0][2]+cur_point[2]],2)
            point_list.pop(0)
            point_list.append(ref_point)
            class_list[0].extend(class_list[1])
            class_list.pop(1)
            break
        # print("distance between 2 cluster center: "+str(cal_dist(point_list[0],cur_point)))
        point_list.append(cur_point)
        print(cur_point)
        i += 1
        if len(class_list)==1:
            len_flag=True

    print('line 196: 5 angstrom aa res seq')
    for i in range(len(point_list)):
        flex_sites=set()
        res={}
        rmsd,pocket=[],[]
        # 5 angstrom aa res seq
        for ligand in class_list[i]:
            print('line 199: ligand->protein: '+str(ligand[0])+'\t'+str(ligand[1][0][0:5]))
            for macro_aa in big_molecules[ligand[1][0][0:5]][0].get_complete_aa():
                if dist_cutoff_cal(ligand[1][1].get_complete_atom(),macro_aa.get_complete_atom(),5)>0:
                    flex_sites.add(macro_aa.get_molecule_resSeq()) 
        
        flex_sites=sorted(list(flex_sites))
        # rmsd cal
        for id in flex_sites:
            # TODO  我认为这里的res_flag没必要设置
            res_flag=False
            max_rmsd=0.0
            if cur_big_molecule.get_complete_aa_by_resSeq(id) == None or cur_big_molecule.get_complete_aa_by_resSeq(id).get_molecule_resName() in ['ALA','GLY','PRO']:
                continue
            for chain_rec in big_molecules:
                template_aa=cur_big_molecule.get_complete_aa_by_resSeq(id)
                compare_aa=big_molecules[chain_rec][0].get_complete_aa_by_resSeq(id)
                if compare_aa is None or template_aa.get_molecule_len() != compare_aa.get_molecule_len() or template_aa.get_molecule_resName()!=compare_aa.get_molecule_resName():
                    continue
                cur_rmsd,atom_list,sc_list,bb_list=compare_aa.rmsd_cal(template_aa)
                res_flag=True
                if cur_rmsd>max_rmsd:
                    max_rmsd=cur_rmsd
            if res_flag:        
                pocket.append(cur_big_molecule.get_complete_aa_by_resSeq(id))
                rmsd.append(max_rmsd)

        # sc 2 lig ctr dist angle cal
        print('line 224: sc 2 lig ctr dist angle cal' )
        for j in range(len(pocket)):
            angle_flag=False
            cur_sc_ctr = geo_ctr(pocket[j],'sc')
            cur_bb_ctr = geo_ctr(pocket[j],'bb')
            cur_dist = cal_dist(cur_sc_ctr,point_list[i])
            cur_angle = 57*angle([point_list[i][0]-cur_bb_ctr[0],point_list[i][1]-cur_bb_ctr[1],point_list[i][2]-cur_bb_ctr[2]],[cur_sc_ctr[0]-cur_bb_ctr[0],cur_sc_ctr[1]-cur_bb_ctr[1],cur_sc_ctr[2]-cur_bb_ctr[2]],True)
            if pocket[i].get_molecule_resName() in ['ARG','LYS','HIS','ASP','GLU','ASN','GLN','LEU','MET','ILE'] and cur_angle<=100 or 360>=cur_angle>=260:
                angle_flag=True      
            elif 90>=cur_angle >0 or 360>=cur_angle>=270:
                angle_flag=True

            res[pocket[j].get_molecule_resSeq()]=[round(rmsd[j],3),round(cur_angle,3),angle_flag,cur_dist]
        
        # remove aa far from lig ctr
        print('line 239: remove aa far from lig ctr' )
        len_res = sorted(res.items(),key=lambda item:item[1][3],reverse=True)
        m=0
        for item in len_res :
            if m >= int(0.2*len(pocket)):
                break
            del res[item[0]]
            m+=1
        res = sorted(res.items(),key=lambda item:item[1][0],reverse=True)
        # output
        prefix = pr_path+'/'+pr+'/flex_site/'+pr
        with open(prefix+'_flexs_'+str(i)+'.txt', 'w') as flex_file, open(prefix+'_pocket_'+str(i)+'.pdb', 'w') as pocket_file:
            n=0
            for item in res:
                flex_file.write(str(item[0])+'\t'+str(item[1][0])+'\t'+str(item[1][1])+'\t'+str(item[1][2])+'\t'+str(item[1][3]))
                if n<5 and item[1][2] == True:
                    flex_file.write('\tcandidate!')
                    n+=1
                flex_file.write('\n')
                
            for aa in pocket:
                pocket_file.writelines(str(aa)+'\n')
            
