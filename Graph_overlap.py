import networkx as nx
import igraph as ig
from math import log
import onmi
import Omega
import time


def my_graph(file_name):
    file = open(file_name)
    g_nx = nx.Graph()
    edge_list = []
    tag = 0
    ind = 0
    for line in file:
        if ind == 0:
            ind += 1
            continue
        line = line.strip('\n')
        num_str = line.split('\t')
        num = [0, 0]
        num[0] = int(num_str[0])
        num[1] = int(num_str[1])
        edge_list.append((num[0] - 1, num[1] - 1))
        if num[1] > tag + 1:
            for I in range(tag + 1, num[1]):
                g_nx.add_node(I - 1)
            tag = num[1]
        g_nx.add_edge(num[0] - 1, num[1] - 1)
    g_ig = ig.Graph(n=tag, edges=edge_list, directed=False)
    return g_nx, g_ig


def graph_comm_dict(community_filename):
    file = open(community_filename)
    comm_dict = {}
    for line in file:
        line = line.strip('\n')
        num_str = line.split(' ')
        num = [0, 0]
        num[0] = int(num_str[0]) - 1
        comms = num_str[1].strip(" ").split(" ")
        comm_dict[num[0]] = [int(i) for i in comms]
    return comm_dict


def graph_comm_list_array(community_file):
    list_array = []
    with open(community_file, 'r') as file:
        for line in file:
            line = line.strip()
            nums = line.split()
            for i in nums:
                j = int(i)
                while len(list_array) < j:
                    list_array.append([])
                list_array[j - 1].append(i)
    return list_array


def graph_comm_set_array(community_filename):
    comm_dict = graph_comm_dict(community_filename)
    keys = comm_dict.keys()
    # print(len(keys))
    length = len(keys)
    set_array = [set() for i in range(length)]
    for i in keys:
        for j in comm_dict[i]:
            set_array[j - 1].add(i)
    res_set_array = [i for i in set_array if i]
    return res_set_array


def graph_real_comm_set(community_filename):
    file = open(community_filename)
    comm_set = []
    for line in file:
        set_list = set()
        line = line.strip('\n')
        num_str = line.split(' ')
        # print(num_str)
        for i in range(len(num_str)):
            set_list.add(int(num_str[i]))
        comm_set.append(set_list)
    return comm_set


def graph_real_comm_list(community_filename):
    file = open(community_filename)
    comm_list = []
    for line in file:
        set_list = list()
        line = line.strip('\n')
        num_str = line.split(' ')
        # print(num_str)
        for i in range(len(num_str)):
            set_list.append(int(num_str[i]))
        comm_list.append(set_list)
    return comm_list


def change_to_mutualdict(A):
    res = dict()
    for index in A.keys():
        if A[index] not in res.keys():
            res[A[index]] = set()
        res[A[index]].add(index)
    return res


def onmi_value(community_comm, community_file):
    # real_comm=graph_comm_set_array(community_file)
    real_comm = graph_real_comm_set(community_file)  # amazon,dblp,youtube时用
    print("数据转化完毕，正在计算onmi值中...")
    t1 = time.time()
    ores = onmi.onmi(community_comm, real_comm)
    t2 = time.time()
    print("计算onmi花费：" + str(t2 - t1), "s")
    return ores


def omega_value(community_comm_omega, community_file):
    real_comm_list = graph_comm_list_array(community_file)
    # real_comm_list = graph_real_comm_list(community_file)
    real_comm_dict_omega = dict()
    for i in range(len(real_comm_list)):
        real_comm_dict_omega[i + 1] = real_comm_list[i]
    print("数据转化完毕，正在计算omega值中...")
    t1 = time.time()
    o = Omega.Omega(community_comm_omega, real_comm_dict_omega).omega_score
    t2 = time.time()
    print("计算omega花费：" + str(t2 - t1), "s")

    return o


def modularity_ov(data, k, res_listoflist):
    t1 = time.time()
    print("数据转化完毕，正在计算重叠模块度中...")
    m = len(k)
    O = [0 for _ in range(m)]
    for i in res_listoflist:  # 计算O的值
        for j in i:
            O[j] += 1
    res = 0.0
    for c in res_listoflist:
        for i in range(len(c)):
            for j in range(i + 1, len(c)):
                temp = data[i][j] - (k[i] * k[j]) / (2 * m)
                temp = temp / (O[i] * O[j])
                res += temp
    t2 = time.time()
    print("重叠模块度计算完毕，花费:", t2 - t1, "s")
    return res / (2 * m)
