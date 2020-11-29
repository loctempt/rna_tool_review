import os
import math


class Atom:
    """
    用于表示PDB文件中的一行（原子）：
    负责实现字符串切片，从表示原子的字符串中读取特定属性。
    提供原子坐标getter以及原子间距离计算函数
    """
    atom_line: str

    def __init__(self, pdb_line: str):
        # 去掉\n
        self.atom_line = pdb_line.strip()

    def __str__(self):
        # 用函数封装对象
        return self.atom_line

    def get_atom(self) -> str:
        if self.atom_line[0] == 'H':
            return str(self.atom_line[0:6])
        else:
            return str(self.atom_line[0:4])

    def get_serial(self) -> int:
        return int(self.atom_line[6:11])

    def get_name(self) -> str:
        return str(self.atom_line[12:16])

    def get_altLoc(self) -> bool:
        return self.atom_line[16]

    def get_resName(self) -> str:
        return str(self.atom_line[17:20])

    def get_chainID(self) -> str:
        return str(self.atom_line[21])

    def get_resSeq(self) -> int:
        return int(self.atom_line[22:26])

    def get_x(self) -> float:
        return float(self.atom_line[30:38])

    def get_y(self) -> float:
        return float(self.atom_line[38:46])

    def get_z(self) -> float:
        return float(self.atom_line[46:54])

    def get_element(self) -> str:
        return str(self.atom_line[76:78])

    def get_temprature_factor(self) -> float:
        return float(self.atom_line[60:66])

    def atom_dist_2(self, another) -> float:
        return (self.get_x() - another.get_x()) ** 2 +\
            (self.get_y() - another.get_y()) ** 2 +\
            (self.get_z() - another.get_z()) ** 2

    def element_compare(self, element: str) -> bool:
        return self.get_element().strip() != element

    def set_resSeq(self, new_resSeq):
        self.atom_line = self.atom_line[:22] + \
            str(new_resSeq).rjust(4, ' ') + self.atom_line[26:]


class AAcid:
    '''
    氨基酸
    负责计算距离、rmsd、删除单聚体多坐标体系
    '''
    atom_list: list

    def __init__(self, chain_lines: list):
        self.atom_list = []
        for line in chain_lines:
            self.atom_list.append(Atom(line))

    def __str__(self):
        return '\n'.join(list(map(lambda atom: str(atom), self.atom_list)))

    def get_complete_atom(self):
        return self.atom_list

    def get_molecule_resName(self):
        return self.atom_list[0].get_resName()

    def get_molecule_resSeq(self):
        return self.atom_list[0].get_resSeq()

    def get_molecule_len(self):
        return len(self.atom_list)

    def min_dist_cal(self, another) -> float:
        '''
        参数：两个分子对象，计算距离
        '''
        min_dist = float('inf')
        for atom in self.atom_list:
            if not atom.element_compare('H'):
                continue
            for another_atom in another.atom_list:
                if not another_atom.element_compare('H'):
                    continue
                cur_dist = math.sqrt(atom.atom_dist_2(another_atom))
                if cur_dist < min_dist:
                    min_dist = cur_dist
        return min_dist

    def rmsd_cal(self, template) -> list:
        '''
        计算氨基酸rmsd，按照pdb记录顺序返回各原子rmsd值
        self是比对链的氨基酸，template是模板链的氨基酸 
        '''
        aminoacid_rmsd = 0.0
        rmsd_atom_list = []
        sidechain_list = []
        backbone_list = []
        for i in range(len(self.atom_list)):
            if not self.atom_list[i].element_compare('H'):
                continue
            cur_val = self.atom_list[i].atom_dist_2(template.atom_list[i])
            if self.atom_list[i].get_name() in ['C', 'N', 'CA', 'O']:
                backbone_list.append(math.sqrt(cur_val))
            else:
                sidechain_list.append(math.sqrt(cur_val))
            rmsd_atom_list.append(math.sqrt(cur_val))
            aminoacid_rmsd += cur_val

        return math.sqrt(aminoacid_rmsd/len(self.atom_list)), rmsd_atom_list, sidechain_list, backbone_list

    def monomer_identify(self):
        cur_list = []
        standard_alt = None
        # 标志位 standard_alt为空时 是false
        for atom in self.atom_list:
            # 记录standard_alt为空的元素
            if atom.get_altLoc() == ' ':
                cur_list.append(atom)
            if standard_alt is None and atom.get_altLoc() != ' ':
                standard_alt = atom.get_altLoc()
                cur_list.append(atom)
            # 若出现第一个标号不为空的：
            elif atom.get_altLoc() == standard_alt:
                cur_list.append(atom)

        if len(cur_list) != len(self.atom_list):
            self.atom_list = cur_list

    def set_resSeq(self, new_resSeq):
        '''
        在序列归一化后重新设置resSeq 
        '''
        for atom in self.atom_list:
            atom.set_resSeq(new_resSeq)


class Chain:
    '''
    用于表示多个氨基酸分子组成的链
    负责生成小分子、区分链号
    '''
    aa_List: list
    aa_fasta_list: list
    chainID: str

    def __init__(self, molecule: list):
        self.aa_List = []
        self.aa_fasta_list = []
        one_aa = []
        cur_resseq = molecule[0][22:26]
        self.chainID = molecule[0][21]
        for line in molecule:
            if len(line) == 0:
                continue
            res_seq = line[22:26]
            if cur_resseq != res_seq:
                self.aa_List.append(AAcid(one_aa))
                one_aa = []
                cur_resseq = res_seq
            one_aa.append(line)

        # 将最后一个分子转换成氨基酸类
        if len(one_aa) > 0:
            self.aa_List.append(AAcid(one_aa))

    def get_chainID(self):
        return self.chainID

    def __str__(self):
        return '\n'.join(list(map(lambda aa: str(aa), self.aa_List)))

    def get_complete_aa_by_resSeq(self, resSeq):
        for aa in self.aa_List:
            if aa.get_molecule_resSeq() == resSeq:
                return aa
        return None

    def get_complete_aa(self):
        return self.aa_List

    def get_complete_atom(self):
        ret = []
        for aa in self.get_complete_aa():
            ret += aa.get_complete_atom()
        return ret

    def pdb_2_fasta(self):
        """
        pdb链转换为fasta
        """
        if len(self.aa_fasta_list) == len(self.aa_List):
            return
        aa_codes = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E','PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N','PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'VAL': 'V', 'TYR': 'Y', 'TRP': 'W', 'THR': 'T','2AS': 'D', '3AH': 'H', '5HP': 'E', 'ACL': 'R', 
        'AIB': 'A', 'ALM': 'A', 'ALO': 'T', 'ALY': 'K','ARM': 'R', 'ASA': 'D', 'ASB': 'D', 'ASK': 'D',
        'ASL': 'D', 'ASQ': 'D', 'AYA': 'A', 'BCS': 'C','BHD': 'D', 'BMT': 'T', 'BNN': 'A', 'BUC': 'C', 
        'BUG': 'L', 'C5C': 'C', 'C6C': 'C', 'CCS': 'C','CEA': 'C', 'CHG': 'A', 'CLE': 'L', 'CME': 'C',
        'CSD': 'A', 'CSO': 'C', 'CSP': 'C', 'CSS': 'C','CSW': 'C', 'CXM': 'M', 'CY1': 'C', 'CY3': 'C', 
        'CYG': 'C', 'CYM': 'C', 'CYQ': 'C', 'DAH': 'F','DAL': 'A', 'DAR': 'R', 'DAS': 'D', 'DCY': 'C', 
        'DGL': 'E', 'DGN': 'Q', 'DHA': 'A', 'DHI': 'H','DIL': 'I', 'DIV': 'V', 'DLE': 'L', 'DLY': 'K', 
        'DNP': 'A', 'DPN': 'F', 'DPR': 'P', 'DSN': 'S','DSP': 'D', 'DTH': 'T', 'DTR': 'W', 'DTY': 'Y', 
        'DVA': 'V', 'EFC': 'C', 'FLA': 'A', 'FME': 'M','GGL': 'E', 'GLZ': 'G', 'GMA': 'E', 'GSC': 'G',
        'HAC': 'A', 'HAR': 'R', 'HIC': 'H', 'HIP': 'H','HMR': 'R', 'HPQ': 'F', 'HTR': 'W', 'HYP': 'P', 
        'IIL': 'I', 'IYR': 'Y', 'KCX': 'K', 'LLP': 'K','LLY': 'K', 'LTR': 'W', 'LYM': 'K', 'LYZ': 'K',
        'MAA': 'A', 'MEN': 'N', 'MHS': 'H', 'MIS': 'S','MLE': 'L', 'MPQ': 'G', 'MSA': 'G', 'MSE': 'M', 
        'MVA': 'V', 'NEM': 'H', 'NEP': 'H', 'NLE': 'L','NLN': 'L', 'NLP': 'L', 'NMC': 'G', 'OAS': 'S', 
        'OCS': 'C', 'OMT': 'M', 'PAQ': 'Y', 'PCA': 'E','PEC': 'C', 'PHI': 'F', 'PHL': 'F', 'PR3': 'C',
        'PRR': 'A', 'PTR': 'Y', 'SAC': 'S', 'SAR': 'G','SCH': 'C', 'SCS': 'C', 'SCY': 'C', 'SEL': 'S', 
        'SEP': 'S', 'SET': 'S', 'SHC': 'C', 'SHR': 'K','SOC': 'C', 'STY': 'Y', 'SVA': 'S', 'TIH': 'A', 
        'TPL': 'W', 'TPO': 'T', 'TPQ': 'A', 'TRG': 'K','TRO': 'W', 'TYB': 'Y', 'TYQ': 'Y', 'TYS': 'Y', 
        'TYY': 'Y', 'AGM': 'R', 'GL3': 'G', 'SMC': 'C', 'ASX': 'B', 'CGU': 'E', 'CSX': 'C', 'GLX': 'Z'
        }
        for cur_aa in self.aa_List:
                self.aa_fasta_list.append(
                    aa_codes[cur_aa.get_molecule_resName()])

    def sequence_match(self, another):
        '''
        self：比对链 another：模板链
        needlman wunsch
        '''
        # TODO: 待查，双人查
        p = 1
        q = 1
        # 初始化 列：mthseq 行：tpseq
        resMth = ""
        resTp = ""
        dp_match = [[0 for i in range(len(self.aa_fasta_list)+1)]
                    for i in range(len(another.aa_fasta_list)+1)]
        # 初始化
        for i in range(len(self.aa_fasta_list)+1):
            dp_match[0][i] = -i
        for i in range(len(another.aa_fasta_list)+1):
            dp_match[i][0] = -i
        # dp_matrix
        while p < len(another.aa_fasta_list)+1:
            q = 1
            while q < len(self.aa_fasta_list)+1:
                if another.aa_fasta_list[p-1] == self.aa_fasta_list[q-1]:
                    dp_match[p][q] = max(
                        dp_match[p-1][q-1]+1, dp_match[p-1][q]-1, dp_match[p][q-1]-1)
                else:
                    dp_match[p][q] = max(
                        dp_match[p-1][q-1]-1, dp_match[p-1][q]-1, dp_match[p][q-1]-1)
                q += 1
            p += 1

        p -= 1
        q -= 1

        # 回溯dp_match
        while p > 0 and q > 0:
            top = dp_match[p - 1][q] - 1
            left = dp_match[p][q - 1] - 1
            if another.aa_fasta_list[p - 1] == self.aa_fasta_list[q - 1]:
                left_top = dp_match[p - 1][q - 1] + 1
                if left_top >= left:
                    if left_top >= top:  # 左上角值最大
                        resMth += self.aa_fasta_list[q-1]
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
                        q -= 1
                    else:  # 上边值最大
                        resMth += "-"
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
                else:
                    if left >= top:  # 左边值最大
                        resMth += self.aa_fasta_list[q-1]
                        resTp += "-"
                        q -= 1
                    else:  # 上边值最大
                        resMth += "-"
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
            elif another.aa_fasta_list[p - 1] != self.aa_fasta_list[q - 1]:
                left_top = dp_match[p - 1][q - 1] - 1
                if left_top >= left:
                    if left_top >= top:  # 左上角值最大
                        resMth += self.aa_fasta_list[q-1]
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
                        q -= 1
                    else:  # 上边值最大
                        resMth += "-"
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
                else:
                    if left >= top:  # 左边值最大
                        resMth += self.aa_fasta_list[q-1]
                        resTp += "-"
                        q -= 1
                    else:  # 上边值最大
                        resMth += "-"
                        resTp += another.aa_fasta_list[p-1]
                        p -= 1
            # print('p:',p,'q:',q)
        while p > 0:
            resMth += "-"
            resTp += another.aa_fasta_list[p-1]
            p -= 1
        while q > 0:
            resTp += "-"
            resMth += self.aa_fasta_list[q-1]
            q -= 1
        return resMth[::-1], resTp[::-1]

    def res_seq_edit(self, mth, tp):
        '''
        对已经序列匹配好的序列 进行编号         '''
        mthList = []
        tpList = []
        tpltIdx = 1
        for idx in tp:
            # print("我是第一个循环：",m)
            if idx == '-':
                tpList.append(0)
            else:
                tpList.append(tpltIdx)
                tpltIdx += 1

        # 找到第一个两两配对的位置
        first_pairing_position = 0
        while first_pairing_position < len(mth):
            if mth[first_pairing_position] == tp[first_pairing_position]:
                break
            first_pairing_position += 1

        for forwardPointer in range(first_pairing_position, len(mth)):
            mthList.append(tpList[forwardPointer])
        for backwardPointer in range(first_pairing_position-1, -1, -1):
            mthList.insert(
                0, backwardPointer-first_pairing_position+tpList[first_pairing_position])

        # 更新chain中各氨基酸的resseq
        amino_acid_cnt = 0
        for i in range(len(mthList)):
            if mth[i] == '-':
                continue
            self.aa_List[amino_acid_cnt].set_resSeq(mthList[i])
            amino_acid_cnt += 1

        return mthList

    def mutation_cal(self, tp, mth, mthList):
        mut_list = []
        print('mth:', len(mth), 'tp:', len(tp), 'mthList:', len(mthList))
        # CHECK len(tp) <= len(mth) ?
        for i in range(len(mth)):
            if mth[i] == '-' or tp[i] == '-':
                continue
            if mth[i] != tp[i]:
                mut_list.append(tp[i]+' '+mth[i]+' '+str(mthList[i])+'\n')
        return mut_list


class Molecule:
    """
    用于表示多个PDB行代表的一个小分子（配体、氨基酸）：
    包含原子对象，提供分子间距计算函数
    """
    # chain_list: list
    chain_dict: dict

    def __init__(self, tmp_chains: list):
        """
        接收atom行，构造Atom对象
        """
        self.chain_dict = {}
        for tmp_chain in tmp_chains:
            # self.chain_list.append(Chain(tmp_chain))
            chain = Chain(tmp_chain)
            self.chain_dict[chain.chainID] = chain

    def __str__(self):
        """
        还原chain_list中的每个Atom对象，得到一个字符串列表，
        再将其拼接成一个完整字符串，返回
        """
        return '\n'.join(list(map(lambda chain: str(chain), list(self.chain_dict.values()))))

    def get_chain(self, id: str):
        return self.chain_dict.get(id)

    def get_complete_chain(self):
        # return self.chain_list
        return list(self.chain_dict.values())

    def get_complete_aa(self):
        ret = []
        for chain in self.get_complete_chain():
            ret += chain.get_complete_aa()
        return ret

    def get_complete_atom(self):
        ret = []
        for aa in self.get_complete_aa():
            ret += aa.get_complete_atom()
        return ret


#       基类构造方法即可
class MacroMolecule(Molecule):
    '''
    作为大分子继承分子类 里面是chain的列表
    '''

    def __init__(self, tmp_chains: list):
        super().__init__(tmp_chains)


class MicroMolecule(Molecule):
    '''
    继承分子类 是活性小分子
    '''

    def __init__(self, tmp_chains: list):
        # 由于pdb传进来的数据没有按照链区分开 在调用super方法之前 需要添加这部分代码
        # print('micro mol tmpChain[0][21]:',tmp_chains[0][21])
        # print('micro mol tmpChain[0]:',tmp_chains[0])

        cur_chainID = tmp_chains[0][21]
        one_chain = []
        chains = []
        for pdb_line in tmp_chains:
            if len(pdb_line) == 0 or pdb_line[17:20] == 'HOH':
                continue
            chainID = pdb_line[21]
            if chainID != cur_chainID:
                cur_chainID = chainID
                chains.append(one_chain)
                one_chain = []
            one_chain.append(pdb_line)

        if len(one_chain) != 0:
            chains.append(one_chain)

        super().__init__(chains)


class PDB:
    macro_molecule: MacroMolecule
    micro_molecule: MicroMolecule

    def __init__(self, lines):
        tmp_molecules = []
        tmp_molecule = []
        status_flag = False
        self.micro_molecule = None
        for line in lines:
            if len(line) == 0:
                continue
            line = line.strip()
            category = line[0:6].strip()
            if not status_flag:
                if category not in ['ATOM', 'HETATM']:
                    continue
                else:
                    tmp_molecule.append(line)
                    status_flag = True
            else:
                if category not in ['ATOM', 'HETATM', 'TER', 'ANISOU']:
                    if tmp_molecule != []:
                        tmp_molecules.append(tmp_molecule)
                        tmp_molecule = []
                    break
                if category == 'TER':
                    tmp_molecules.append(tmp_molecule)
                    tmp_molecule = []
                elif category == 'ANISOU':
                    continue
                else:
                    tmp_molecule.append(line)

        if tmp_molecule != []:
            tmp_molecules.append(tmp_molecule)
        flag = False
        # print(tmp_molecules[-1][0][0])
        for mol in tmp_molecules[-1]:
            if mol[0][0] == 'A':
                flag = True
                break
        if flag:
            self.macro_molecule = MacroMolecule(tmp_molecules)
            # print('has macro_molecule')
        else:
            # print('has micro_molecule,macro_molecule')
            # print(tmp_molecules)
            if len(tmp_molecules) == 2:
                self.macro_molecule = MacroMolecule([tmp_molecules[0]])
                self.micro_molecule = MicroMolecule(tmp_molecules[1])
            else:
                macro_mol = []
                for i in range(len(tmp_molecules)-1):
                    macro_mol.append(tmp_molecules[i])
                self.macro_molecule = MacroMolecule(macro_mol)
                self.micro_molecule = MicroMolecule(tmp_molecules[-1])

    def __str__(self):
        return '\n'.join(list(map(
            lambda mol: str(mol),
            self.macro_molecule.get_complete_chain()+self.micro_molecule.get_complete_chain()
        )))

    def molecule_cat(self, id: str):
        if self.micro_molecule != None and self.micro_molecule.get_chain(id) != None:
            # print(str(self.micro_molecule.get_chain(id)))
            return (self.macro_molecule.get_chain(id), self.micro_molecule.get_chain(id))
        else:
            return (self.macro_molecule.get_chain(id), None)

    def get_micro_molecule(self):
        if self.micro_molecule != None:
            return self.micro_molecule
        else:
            return None

    def get_macro_molecule(self):
        return self.macro_molecule
