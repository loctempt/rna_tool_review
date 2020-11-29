from rna_tool_new import *
import os
import shutil
import subprocess
import math
import math
import util
from config import Config


def dist_cutoff_cal(micro_atom_list, macro_atom_list, cutoff):
    """
    计算micro_atom_list链条中，与macro_atom_list链条最小距离不大于cutoff的micro_atom总数
    """
    atom_cnt = 0
    for micro_atom in micro_atom_list:
        if not micro_atom.element_compare('H'):
            continue
        for macro_atom in macro_atom_list:
            if not macro_atom.element_compare('H'):
                continue
            if math.sqrt(micro_atom.atom_dist_2(macro_atom)) > cutoff:
                continue
            atom_cnt += 1
            break
    return atom_cnt


def input_legalization(marked_words, legal_num, case_status):
    flag = True
    while flag is True:
        tmp = input(marked_words)
        if len(tmp) == legal_num and case_status(tmp):
            flag = False
    return tmp


def call_perl_and_rename(rename_to, matchPath, script_path):
    print('start running perl script...')
    os.chdir(matchPath)
    print('call_perl:', os.getcwd())
    ret = subprocess.call(['perl', script_path+'/superimopose.pl'])
    if ret == 0:
        os.rename(matchPath+'/align_tmp.pdb', matchPath+'/'+rename_to)
        print('success')
    else:
        print('failed')


def call_superimopose(pdb_id, data_pth):
    pdb_list = os.listdir(data_pth+'/pre_dealing')
    shutil.copyfile(data_pth+'/pre_dealing/'+pdb_id+'.pdb',
                    data_pth+'/superimpose/template.pdb')
    for i in pdb_list:
        shutil.copyfile(data_pth+'/pdb/'+i, data_pth+'/superimpose/align.pdb')
        call_perl_and_rename(i, data_pth+'/superimpose',
                             script_pth)  # FIXME script_pth未定义
        # TODO script path undefined
        os.remove(data_pth+'/superimpose/align.pdb')
    os.remove(data_pth+'/superimpose/template.pdb')


def prepare_fasta(template_id, pdb: PDB, cluster_dir):
    """
    将pdb中的大分子转成fasta格式，写入缓存文件
    """
    macro_molecule = pdb.get_macro_molecule()
    chains = macro_molecule.get_complete_chain()
    out_file = os.path.join(cluster_dir, template_id+".fasta")
    with open(out_file, 'w'):
        for chain in chains:
            chain.pdb_2_fasta()
        out_file.write('\n'.join(chains))


def get_chain_id_by_cluster(template_id, cluster_dir):
    """
    从聚类结果中找到长度最接近receptor的模板链
    """
    cluster_file = os.path.join(cluster_dir, "firstRes.clstr")
    cluster = util.get_cluster_by_id(
        cluster_file, util.ClsrReader.RECEPTOR_FLAG, util.ClsrReader.RECEPTOR_FLAG)
    length_of_receptor = util.length_of_chain(
        cluster, util.ClsrReader.RECEPTOR_FLAG, util.ClsrReader.RECEPTOR_FLAG)
    cluster = list(filter(
        lambda x: x[1] != util.ClsrReader.RECEPTOR_FLAG and x[1] == template_id, cluster))
    cluster = sorted(cluster, key=lambda x: abs(x[0] - length_of_receptor))
    chain_id = cluster[0][2]
    return chain_id


def clustering_stage_a(template_id, file_dir, cluster_dir):
    """
    模板链ID未知，通过stage_a聚类，找出模板链ID
    """
    obabel_input = os.path.join(file_dir, "dude", "receptor.pdb")
    obabel_output = os.path.join(cluster_dir, "receptor.fasta")
    obabel_cmd = "obabel -ipdb {} -ofasta -O {}".format(
        obabel_input, obabel_output)
    os.system(obabel_cmd)

    cat_input = obabel_output
    cat_output = os.path.join(cluster_dir, template_id+".fasta")
    os.system("cat {} >> {}".format(cat_input, cat_output))

    cdhit_input = cat_output
    cdhit_output = os.path.join(cluster_dir, "firstRes")
    os.system("cdhit -i {} -o {} -c 0.70 -n 5".format(cdhit_input, cdhit_output))


def clustering_stage_b(cluster_dir):
    """
    正式聚类过程
    """
    cdhit_input = os.path.join(cluster_dir, "fasta.txt")
    cdhit_output = os.path.join(cluster_dir, "secondRes")
    os.system("cdhit -i {} -o {} -c 0.70 -n 5".format(cdhit_input, cdhit_output))


def cluster_sanitize(template_id, chain_id, cluster_dir):
    """
    去除长度小于模板链长度50%的短链
    """
    cluster_file = os.path.join(cluster_dir, "secondRes.clstr")
    cluster = util.get_cluster_by_id(cluster_file, template_id, chain_id)
    length_of_template = util.length_of_chain(cluster, template_id, chain_id)
    cluster = list(filter(lambda x: x[0] >= 0.5 * length_of_template, cluster))
    return cluster


def clustering(template_id, file_dir):
    """
    聚类功能

    Returns:
    模板链的chain id，聚类结果（排除与模板链长度差值超过50%的链）
    """
    cluster_dir = os.path.join(file_dir, "序列聚类")
    template_pdb = None
    chain_id = None
    # 构造PDB对象
    with open(os.path.join(file_dir, 'data', template_id+".pdb"), 'r') as pdb_file:
        template_pdb = PDB(pdb_file.readlines())
    # 针对单链/多链分别采取对策
    if util.is_single_chain(template_pdb):
        chain_id = util.get_chain_id(template_pdb)
    else:
        prepare_fasta(template_id, template_pdb, cluster_dir)
        clustering_stage_a(template_id, file_dir, cluster_dir)
        chain_id = get_chain_id_by_cluster(template_id, cluster_dir)
    clustering_stage_b(cluster_dir)
    sanitized_cluster = cluster_sanitize(template_id, chain_id, cluster_dir)
    return chain_id, sanitized_cluster


target_dict = util.read_targets(Config.TARGETS_FILE_PATH)


def get_template_id(protein_name):
    return target_dict[protein_name.upper()]


# dirPath = input("输入当前这一类蛋白（如激酶）的目录:")
dirPath = Config.DIR_PATH
dirList = os.listdir(dirPath)

for pr_class_name in dirList:
    # os.system('pause')
    filePath = os.path.join(dirPath, pr_class_name)
    # if not os.path.isdir(filePath) or pr_class_name not in ['braf']:
    #     continue

    fileList = os.listdir(os.path.join(filePath, 'data'))
    # PDB_dict = {}
    templateID = get_template_id(pr_class_name)
    template_chain_id, cluster = clustering(templateID, filePath)

    rough_sele_pr_dict = {}
    for file in fileList:
        print(file)
        #  TODO 将两个部分交换后 加判断 以确定file在不在列表中，
        if file not in list:
            continue
        with open(os.path.join(filePath, 'data', file), 'r') as cur_file:
            cur_PDB = PDB(cur_file.readlines())
            # 将大小分子删去单聚体多坐标体系
            for aa in cur_PDB.macro_molecule.get_complete_aa():
                aa.monomer_identify()
            # 对含有小分子的蛋白质做以下操作
            if cur_PDB.micro_molecule != None:
                for aa in cur_PDB.micro_molecule.get_complete_aa():
                    aa.monomer_identify()

                # 功能： 删去没有3个碳、距离与蛋白质在2埃以内的小分子
                for chain in cur_PDB.micro_molecule.get_complete_chain():
                    tmp_2_chain = []
                    for aa in chain.get_complete_aa():
                        carbon_cnt = 0
                        for atom in aa.get_complete_atom():
                            if not atom.element_compare('C'):
                                carbon_cnt += 1
                        if carbon_cnt < 3:
                            continue
                        min_dist = float('inf')
                        macro_chain = cur_PDB.macro_molecule.get_chain(
                            chain.get_chainID())
                        for macro_aa in macro_chain.get_complete_aa():
                            cur_dist = macro_aa.min_dist_cal(aa)
                            if min_dist > cur_dist:
                                min_dist = cur_dist
                        # print('min_dist',min_dist)
                        if carbon_cnt >= 3 and min_dist > 2.0:
                            tmp_2_chain.append(aa)

                    chain.aa_List = tmp_2_chain
                    # remove导致迭代不会经过下一个元素
                    # 因为遍历会每次取下一个元素
                    # 而remove使下一个元素前移
                    # 引用不传递！！！！！

            macro_chain, micro_chain = cur_PDB.molecule_cat(curChainID)
            if micro_chain is not None:
                for rough_micro_mol in micro_chain.get_complete_aa():
                    one_Mol_total_cutoff = dist_cutoff_cal(
                        rough_micro_mol.get_complete_atom(), macro_chain.get_complete_atom(), 5)
                if float(one_Mol_total_cutoff)/len(rough_micro_mol.get_complete_atom()) < 0.6:
                    micro_chain.aa_List.remove(rough_micro_mol)
                rough_sele_pr_dict[str(
                    curPdb+curChainID)] = [macro_chain, micro_chain]
            else:
                rough_sele_pr_dict[str(curPdb+curChainID)] = [macro_chain]

            # PDB_dict[str(file)[:4]] = cur_PDB

    # =========================
    # 功能：读outRes.clsr文件
    # =========================
    # 输出pdb链
    # rough_sele_pr_dict = {}
    # clsrNum = input("plz input "+pr_class_name+" template class ID:")

    # with open(os.path.join(filePath,'序列聚类','outRes.clstr'), 'r') as clsrFile:
    #     clsrnum_found = False
    #     maxchainlen = -1
    #     for line in clsrFile.readlines():
    #         line = line.split()
    #         if line[0] == '>Cluster' and line[1] == str(clsrNum):
    #             clsrnum_found = True
    #             # 跳过目标段落第一行
    #             continue
    #         if clsrnum_found:
    #             if line[0] == '>Cluster':
    #                 # 来到下一段落，结束循环
    #                 break
    #             curPdb = line[2][1:5]
    #             curChainID = line[2][6]

                # macro_chain, micro_chain = PDB_dict[curPdb.lower()].molecule_cat(
                #     curChainID)
                # if len(macro_chain.get_complete_aa()) >= maxchainlen:
                #     maxchainlen = len(macro_chain.get_complete_aa())
                # if micro_chain is not None:
                #     rough_sele_pr_dict[str(
                #         curPdb+curChainID)] = [macro_chain, micro_chain]
                # else:
                #     rough_sele_pr_dict[str(curPdb+curChainID)] = [macro_chain]

    # =========================
    # 筛选合适的肽链以及小分子
    # =========================
    # 将小于0.55长度的链都删去
    # 并在剩余的rough_dict中筛选合适的小分子 距离小于5 60%
    # tmp_dict = {}
    # for key, value in rough_sele_pr_dict.items():
    #     macro_chain = value[0]
    #     micro_chain = value[1] # if len(value) == 2 else None
    #     if float(len(macro_chain.get_complete_aa()))/maxchainlen > 0.5:
            # if len(value) == 1:
            #     tmp_dict[key] = value
            #     continue
    #         for rough_micro_mol in micro_chain.get_complete_aa():
    #             one_Mol_total_cutoff = dist_cutoff_cal(
    #                 rough_micro_mol.get_complete_atom(), macro_chain.get_complete_atom(), 5)
    #             if float(one_Mol_total_cutoff)/len(rough_micro_mol.get_complete_atom()) < 0.6:
    #                 micro_chain.aa_List.remove(rough_micro_mol)
    #         tmp_dict[key] = value
    # rough_sele_pr_dict = tmp_dict

    # 输出单链文件 检查小分子和肽链 当前筛选的结果是否合适 属于中间结果 可删
    # for key, value in rough_sele_pr_dict.items():
    #     with open(os.path.join(filePath,'output',key+'.pdb'), 'w') as outFile:
    #         outFile.writelines(str(value[0]))
    #         outFile.write('\nTER\n')
    #         if len(value) > 1:
    #             outFile.writelines(str(value[1]))

    # =========================
    # 序列叠合编号 鉴定突变位点
    # =========================
    # TODO template ID 是dude 模板蛋白 id
    template_chain = rough_sele_pr_dict[templateID][0]
    template_chain.pdb_2_fasta()
    for key, value in rough_sele_pr_dict.items():
        macro_chain = value[0]
        macro_chain.pdb_2_fasta()
        mth, tp = macro_chain.sequence_match(template_chain)
        mthList = macro_chain.res_seq_edit(mth, tp)
        # 找突变点, 暂时未用上
        # with open(os.path.join(filePath,'mutations',key+'.txt'), 'w') as mutfile:
        #     for i in macro_chain.mutation_cal(tp, mth, mthList):
        #         mutfile.write(str(i))

    # =========================
    # 输出编号归一的pdb (输出时大小分子间要插入TER)
    # =========================
    for key, value in rough_sele_pr_dict.items():
        macro_chain = value[0]
        micro_chain = value[1] if len(value) == 2 else None
        with open(os.path.join(filePath, 'pre_dealing', key+'.pdb'), 'w') as outFile:
            outFile.writelines(str(macro_chain))
            outFile.write('\nTER\n')
            if len(value) > 1:
                outFile.writelines(str(micro_chain))

    # superimpose
    # print('下一步：蛋白叠合，请移至服务器完成')
    call_superimopose(templateID, dirPath+'/'+pr_class_name)
    os.system('pause')

    # =========================
    # 计算口袋
    # =========================
    # 读文件构造pdb对象
    # 计算与各个小分子距离在5埃内的分子  构造大小分子dict{key,[macro,micro]}
    #  没有小分子的直接输出一个大分子即可 大小分子间加TER
    # superimposefileList = os.listdir(os.path.join(filePath,'superimpose','input'))
    # # rough_sele_pr_dict.clear()
    # rough_sele_pr_dict = {}
    # for surfile in superimposefileList:
    #     print(surfile)
    #     with open(os.path.join(filePath,'superimpose','input',surfile), 'r') as cur_file:
    #         cur_PDB = PDB(cur_file.readlines())
    #         if cur_PDB.micro_molecule != None:
    #             rough_sele_pr_dict[str(surfile[:5])] = [cur_PDB.macro_molecule.get_chain(
    #                 surfile[4]), cur_PDB.micro_molecule.get_chain(surfile[4])]
    #         else:
    #             rough_sele_pr_dict[str(surfile[:5])] = [
    #                 cur_PDB.macro_molecule.get_chain(surfile[4])]

    # specified_sele_pkt_dict = {}
    # ligand_cnt = 0
    # tmp_key = 0
    # for key, value in rough_sele_pr_dict.items():
    #     macro_chain = value[0]
    #     micro_chain = value[1] if len(value) == 2 else None
    #     if len(value) == 1:
    #         specified_sele_pkt_dict[key] = [macro_chain]
    #     else:
    #         # 对每种key单独编号
    #         if tmp_key != key:
    #             tmp_key = key
    #             ligand_cnt = 0
    #         for micro_aa in micro_chain.get_complete_aa():
    #             tmp_pkt_list = []
    #             for macro_aa in macro_chain.get_complete_aa():
    #                 dist_cutoff = dist_cutoff_cal(
    #                     macro_aa.get_complete_atom(), micro_aa.get_complete_atom(), 5)
    #                 if dist_cutoff == 0:
    #                     continue
    #                 # 记录落入口袋半径的氨基酸
    #                 tmp_pkt_list.append(macro_aa)
    #             specified_sele_pkt_dict[key+'_' +
    #                                     str(ligand_cnt)] = [tmp_pkt_list, micro_aa]
    #             ligand_cnt += 1

    # for key, value in specified_sele_pkt_dict.items():
    #     macro_chain = value[0]
    #     micro_chain = value[1] if len(value) == 2 else None
    #     print(key)
    #     with open(os.path.join(filePath,'superimpose','output',key+'.pdb'), 'w')as outFile:
    #         if len(value) == 1:
    #             outFile.writelines(str(macro_chain))
    #         else:
    #             for i in macro_chain:
    #                 outFile.writelines(str(i))
    #                 outFile.write('\n')
    #             outFile.write('TER\n')
    #             outFile.write(str(micro_chain))

    # pymol 删除不可用的口袋 (手动处理后放至文件夹：flex_site-> input)
    # print('请将候选口袋初步判断后，将候选口袋移至flex_site/input文件夹，进行下一步')
    # os.system('pause')
    # # =========================
    # # 鉴定柔性
    # # ========================
    # # 读文件构造pdb对象
    # # set()存各个现存的氨基酸编号  遍历上一步的dict中取通用口袋
    # # 计算rmsd值 每一个氨基酸区分和不区分 sidechain 和 backbone分别对原子计算
    # # 找到sidechain和backbone的最大值 实际上要找到每个氨基酸的最大rmsd值 保存该口袋
    # # 格式：2OF2A
    # rough_template_pkt = input_legalization(
    #     "DUD-E" + pr_class_name + " receptor(2OF2A)：", 5, str.isupper)
    # rough_pkfileList = os.listdir(os.join(filePath,'flex_site','input'))
    # specified_sele_pkt_dict = {}
    # # specified_sele_pkt_dict.clear()
    # for rgh_file in rough_pkfileList:
    #     with open(os.path.join(filePath,'flex_site','input',rgh_file), 'r') as cur_file:
    #         print(rgh_file)
    #         cur_PDB = PDB(cur_file.readlines())
    #         if cur_PDB.micro_molecule != None:
    #             print(rgh_file[4])
    #             # print(str(cur_PDB.macro_molecule.get_chain(rgh_file[4])))
    #             # if not specified_sele_pkt_dict.get(str(rgh_file[:5])):
    #             #     specified_sele_pkt_dict[str(rgh_file[:5])]=cur_PDB.macro_molecule.get_chain(rgh_file[4])
    #             # else:
    #             specified_sele_pkt_dict[str(
    #                 rgh_file[:5])] = cur_PDB.macro_molecule.get_chain(rgh_file[4])

    # pkt_resSeq = set()
    # for key, value in specified_sele_pkt_dict.items():
    #     for macro_aa in value.get_complete_aa():
    #         pkt_resSeq.add(macro_aa.atom_list[0].get_resSeq())

    # pkt_resSeq = sorted(list(pkt_resSeq))

    # specified_sele_pkt_dict.clear()
    # for key, value in rough_sele_pr_dict.items():
    #     tmp_sp_pkt_list = []
    #     for macro_aa in value[0].get_complete_aa():
    #         if macro_aa.atom_list[0].get_resSeq() in pkt_resSeq:
    #             tmp_sp_pkt_list.append(macro_aa)
    #     specified_sele_pkt_dict[key] = [tmp_sp_pkt_list]

    # template_pkt = specified_sele_pkt_dict[rough_template_pkt][0]
    # max_rmsd_pkt = dict(map(lambda item: (item, float(0)), pkt_resSeq))
    # max_pkt = dict(map(lambda item: (item, None), pkt_resSeq))
    # for key, value in specified_sele_pkt_dict.items():
    #     print(key, 'len:', len(value[0]))
    #     for cur_aa_no in range(len(value[0])):
    #         cur_resSeq = value[0][cur_aa_no].get_molecule_resSeq()
    #         cur_resName = value[0][cur_aa_no].get_molecule_resName()
    #         cur_aa_len = value[0][cur_aa_no].get_molecule_len()
    #         for temp_aa in template_pkt:
    #             if temp_aa.get_molecule_resSeq() != cur_resSeq or temp_aa.get_molecule_resName() != cur_resName or len(temp_aa.atom_list) != cur_aa_len:
    #                 continue
    #             rmsd, atom_list, sc_list, bb_list = value[0][cur_aa_no].rmsd_cal(
    #                 temp_aa)

    #             if rmsd >= max_rmsd_pkt[value[0][cur_aa_no].atom_list[0].get_resSeq()]:
    #                 max_rmsd_pkt[cur_resSeq] = rmsd
    #                 max_pkt[cur_resSeq] = temp_aa
    # # 删除max_pkt中为None的aa记录
    # for key in list(max_pkt.keys()):
    #     if max_pkt[key] is None:
    #         del max_pkt[key]
    #         del max_rmsd_pkt[key]

    # # 添加计算中心点距离的功能

    # # 添加rmsd值排序以及标注氨基酸的功能
    # sorted_max_rmsd_pkt = sorted(
    #     max_rmsd_pkt.items(), key=lambda item: item[1], reverse=True)

    # with open(os.path.join(filePath,'flex_site','output',pr_class_name+'_'+rough_template_pkt+'_pocket.pdb'), 'w') as pocket, open(os.path.join(filePath,'flex_site','output',rough_template_pkt+'_pocket_rmsd.txt'), 'w') as pocket_rmsd:
    #     # for i in range(len,ax)
    #     for pkt_key in max_pkt:
    #         pocket.writelines(str(max_pkt[pkt_key]))
    #         pocket.write('\n')
    #     for rmsd_key in sorted_max_rmsd_pkt:
    #         pocket_rmsd.write(str(rmsd_key[0])+'  '+str(rmsd_key[1]))
    #         # print(max_pkt[rmsd_key[0]].get_molecule_resName())
    #         if max_pkt[rmsd_key[0]].get_molecule_resName() in ['ALA', 'GLY', 'PRO']:
    #             pocket_rmsd.write(
    #                 '  '+max_pkt[rmsd_key[0]].get_molecule_resName())
    #         pocket_rmsd.write('\n')
    #         # print(max_rmsd_pkt[key],'\n',max_pkt[key])
