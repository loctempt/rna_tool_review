import os
import openpyxl.writer.excel
from matplotlib import pyplot as plt 
import matplotlib
import numpy as np
import math

def cal_y(x,numlist,lig_num):
    y=np.arange(0,100,0.001)
    for i in range(len(x)):
        sub_arr = numlist[:int(x[i] / 100 * len(numlist))]
        ligand_cnt=0
        for j in sub_arr:
            if j=='ligand':
                ligand_cnt+=1
        if len(sub_arr)!=0:
            y[i]=(float(ligand_cnt)/lig_num)*100
    return y


inputPath=r'C:\Users\lenovo\Desktop\biyesheji'
protein_type,protein=map(str,input('please input protein type & name: \n').split())
workbook = openpyxl.Workbook()

fig = plt.figure(num=1, figsize=(15, 8),dpi=80) 
plt.title(protein+" Roc Curve") 
plt.xlabel("Decoy Founds%")
plt.ylabel("Ligand Founds%")
lig_cnt=0

# itero = 0
for iter in ['flex','rigid']:   
    print(iter)
    worksheet = workbook.create_sheet(protein+'_'+iter,0)
    res_list = [] 
    for jter in ['ligand','decoy']:
        print(jter)
        molecule_list = os.listdir(inputPath+'\\'+protein_type+'\\'+protein+'\\dock\\output\\'+str(iter)+'\\'+str(jter))
        if jter=='ligand':
            lig_cnt=len(molecule_list)
        for molecule in molecule_list:
            # print(molecule)
            
            if not os.path.isdir(inputPath+'\\'+protein_type+'\\'+protein+'\\dock\\output\\'+str(iter)+'\\'+str(jter)+'\\'+str(molecule)):
                continue
            with open(inputPath+'\\'+protein_type+'\\'+protein+'\\dock\\output\\'+str(iter)+'\\'+str(jter)+'\\'+str(molecule)+'\\log.txt','r') as cur_file:
                print(str(molecule))
                dock_score = min([float(record.split()[1]) for record in cur_file.readlines()[25:-1]])
                res_list.append([jter,molecule,dock_score])

    res_list = sorted(res_list,key=(lambda x: x[2]))
    x = np.arange(0,100,0.001) 
    y = cal_y(x,[rec[0] for rec in res_list],lig_cnt)
    plt.plot(x,y,label=iter) 

    cnt=1
    for res_record in res_list:
        worksheet.cell(cnt,1).value=res_record[0]
        worksheet.cell(cnt,2).value=res_record[1]
        worksheet.cell(cnt,3).value=res_record[2]
        cnt+=1

workbook.save(inputPath+'\\'+protein_type+'\\'+protein+'\\dock\\output\\'+protein+'.xlsx')
plt.legend()
plt.show()




# decoyresultfilePath=r'C:\Users\lenovo\Desktop\biyesheji\akt1\virtual screen\autodockVina\output\rigid\decoy\result'
# ligandresultfilePath=r'C:\Users\lenovo\Desktop\biyesheji\akt1\virtual screen\autodockVina\output\rigid\ligand\result'
# decoyMolsPath=r'C:\Users\lenovo\Desktop\biyesheji\akt1\virtual screen\autodockVina\dockprep\decoy\mol2'
# ligandMolsPath=r'C:\Users\lenovo\Desktop\biyesheji\akt1\virtual screen\autodockVina\dockprep\ligand\mol2'
# outPath=r'C:\Users\lenovo\Desktop\biyesheji\akt1\virtual screen\autodockVina\output\rigid'

# decoyList=os.listdir(decoyresultfilePath)
# ligandList=os.listdir(ligandresultfilePath)
# workbook = openpyxl.Workbook()
# worksheet = workbook.create_sheet('akt1',0)
# tempdict={}

# for decoy in decoyList:
#     decoy=decoy.split('.')
#     with open(decoyresultfilePath+'\\'+decoy[0]+'.'+decoy[1],'r') as decoyresFile,open(decoyMolsPath+'\\'+decoy[0]+'.mol2','r')as decoymol2File:
#         templist=[]
#         for line in decoyresFile.readlines():
#             line=line.split()
#             templist.append('decoy')
#             templist.append(int(decoy[0]))
#             templist.append(float(line[1]))
            

#         for line in decoymol2File.readlines():
#             if line[0:4]=='ZINC':
#                 if tempdict.__contains__(str(line)):
#                     if tempdict[str(line)][2]>templist[2]:
#                         tempdict[str(line)]=templist
#                 else:
#                     tempdict[str(line)]=templist
#                 break
            

# for ligand in ligandList:
#     ligand=ligand.split('.')
#     with open(ligandresultfilePath+'\\'+ligand[0]+'.result','r') as ligandresFile,open(ligandMolsPath+'\\'+ligand[0]+'.mol2','r') as ligandmol2File:
#         templist=[]
#         for line in ligandresFile.readlines():
#             line=line.split()
#             templist.append('ligand')
#             templist.append(int(ligand[0]))
#             templist.append(float(line[1]))
            

#         for line in ligandmol2File.readlines():
#             if line[0:6]=='CHEMBL':
#                 if tempdict.__contains__(str(line)) :
#                     if tempdict[str(line)][2]>templist[2]:
#                         tempdict[str(line)]=templist
#                 else:
#                     tempdict[str(line)]=templist
#                 break
            


# cnt=1
# print(len(tempdict))

# for key in tempdict:
#     worksheet.cell(cnt,1).value=key
#     temp=tempdict[key]
#     print(key,cnt,tempdict[key][1])
#     worksheet.cell(cnt,2).value=tempdict[key][0]
#     worksheet.cell(cnt,3).value=tempdict[key][1]
#     worksheet.cell(cnt,4).value=tempdict[key][2]
#     cnt+=1
# workbook.save(outPath+'\\akt1_rigid.xlsx')

