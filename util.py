from rna_tool_new import *


def read_targets(target_file='targets.out'):
    """
    读targets文件，形成 “蛋白类别=>模板蛋白” 的映射关系
    """
    target_dict = {}
    with open(target_file) as targets:
        lines = targets.readlines()
        for line in lines:
            line = line.split()
            target_dict[line[0]] = line[1]
    return target_dict


def get_complete_chain(pdb: PDB):
    """
    取出pdb文件中所有大分子链
    """
    macro_molecule = pdb.get_macro_molecule()
    complete_chain = macro_molecule.get_complete_chain()
    return complete_chain


def get_chain_id(pdb: PDB):
    """
    取单链蛋白的链id
    """
    complete_chain = get_complete_chain(pdb)
    chain = complete_chain[0]
    return chain.get_chainID()


def is_single_chain(pdb: PDB):
    """
    检查给定pdb是否为单链结构
    """
    complete_chain = get_complete_chain(pdb)
    return len(complete_chain) == 1


def length_of_chain(cluster: list, pdb_id, chain_id):
    """
    返回给定肽链的长度（aa）
    若给定肽链不存在，返回-1
    pdb_id, chain_id同时传入ClsrReader.RECEPTOR_FLAG，就是找receptor
    """
    for item in cluster:
        if chain_id == ClsrReader.CHAIN_ID_ANY and item[1] == pdb_id.upper():
            return item[0]
        if item[1] == pdb_id.upper() and item[2] == chain_id.upper():
            return item[0]
    else:
        return -1


def get_cluster_by_id(cluster_file, pdb_id, chain_id):
    """
    从clusters中取出包含pdb_id和chain_id的cluster
    pdb_id, chain_id同时传入ClsrReader.RECEPTOR_FLAG，就是找receptor
    """
    with ClsrReader(cluster_file) as clsr_reader:
        while True:
            cluster = clsr_reader.next_cluster()
            if cluster is None:
                return None
            length_of_receptor = length_of_chain(
                cluster, pdb_id.upper(), chain_id.upper())
            if length_of_receptor > -1:
                return cluster


def count_clusters(cluster_file):
    """
    统计cluster_file中包含的cluster数量
    """
    cnt = 0
    with ClsrReader(cluster_file) as clsr_reader:
        while True:
            cluster = clsr_reader.next_cluster()
            if cluster is None:
                break
            cnt += 1
    return cnt


class ClsrReader:
    RECEPTOR_FLAG = "receptor"
    CHAIN_ID_ANY = "any_chain"

    def __init__(self, clsr_file):
        self.clsr_file = open(clsr_file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.clsr_file.close()

    def _parse_cluster(self, revert=0):
        """
        将聚类得到的类别信息转换为list
        """
        # Once a parsing progress is done, the next ">Cluster" line
        # is skipped by finding the ending point of current cluster,
        # and the next "1 100aa, >6PS2:A..." line is skipped by caller
        # while checking EOF sign. So, whenever a line containing
        # chain data is skipped, we need to revert the file pointer
        # to retain this line.
        file_pos = self.clsr_file.tell()
        self.clsr_file.seek(file_pos - revert)
        cluster = []
        while True:
            line = self.clsr_file.readline()
            if not line or line.startswith(">Cluster"):
                # 读完一个类，或者读到文件结尾时返回
                return cluster
            values = line.split()
            # original format: "123aa,".
            chain_length = int(values[1][:-3])
            # 为clustering_stage_a标记模板链
            if values[2].startswith(">receptor"):
                pdb_id = ClsrReader.RECEPTOR_FLAG
                chain_id = ClsrReader.RECEPTOR_FLAG
            # 正常读取
            else:
                pdb_id = values[2][1:5]
                chain_id = values[2][6]
            cluster.append((chain_length, pdb_id, chain_id))

    def next_cluster(self):
        """
        读聚类文件，从中取出下一个类
        """
        while True:
            line = self.clsr_file.readline()
            if not line:
                return None
            if line.startswith(">Cluster"):
                return self._parse_cluster()
            return self._parse_cluster(len(line))
