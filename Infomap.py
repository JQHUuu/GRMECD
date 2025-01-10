import infomap
import collections
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import Graph_overlap
from math import log
import numpy as np

import main_tool_new


def get_reconstructed_graph(resourceMatrix, num_edge=1):
    G1 = nx.Graph()
    for i in range(len(resourceMatrix)):
        max_key = [0] * num_edge
        max_value = [0] * num_edge
        temp_value = resourceMatrix[i][i]
        resourceMatrix[i][i] = 0
        for j in range(num_edge):
            max_key[j] = np.argmax(resourceMatrix[i])
            max_value[j] = resourceMatrix[i][max_key[j]]
            resourceMatrix[i][max_key[j]] = 0
        for j in range(num_edge):
            resourceMatrix[i][max_key[j]] = max_value[j]
            G1.add_edge(i, max_key[j])
        resourceMatrix[i][i] = temp_value
    return G1


class Graph:
    graph = nx.DiGraph()

    def __init__(self):
        self.graph = nx.DiGraph()

    def createGraph(self, filename):
        file = open(filename, 'r')

        for line in file.readlines():
            nodes = line.split()
            edge = (int(nodes[0]), int(nodes[1]))
            self.graph.add_edge(*edge)

        return self.graph


def get_G(filename):
    G_node = set()
    G_edge1 = []
    G_edge2 = []
    with open(filename) as f:
        G = nx.Graph()
        index = 0
        for line in f:
            line = line.split(' ')
            u = int(line[0])
            v = int(line[1])
            u = u  # 从1开始要-1
            v = v
            G_node.add(u)
            G_edge1.append(u)
            G_node.add(v)
            G_edge2.append(v)
    f.close()
    G_node = list(G_node)
    G_node = sorted(G_node)
    for n in range(len(G_node)):
        G.add_node(G_node[n])
    for e in range(len(G_edge1)):
        G.add_edge(G_edge1[e], G_edge2[e])
    return G


def get_G_gml(filename):
    file_read = open(filename, "r")
    real_label1 = []
    nodes = []
    edges1 = []
    edges2 = []
    for line in file_read:
        line = line.split('\t')
        nodes.append(int(line[0]))
        for i in range(2, len(line) - 1):
            edges1.append(int(line[0]))
            edges2.append(int(line[i].split("(")[0]))
        real_label1.append(line[len(line) - 1].split("\n")[0])
    G = nx.DiGraph()
    for n in range(len(nodes)):
        G.add_node(nodes[n])
    for e in range(len(edges1)):
        G.add_edge(edges1[e], edges2[e])
    return G


class Infomap:
    graph = Graph()
    n = len(graph.graph.nodes)

    def __init__(self, G):
        self.graph = G
        self.n = len(G.nodes)

    def findCommunities(self, G):
        infomapWrapper = infomap.Infomap("--two-level --directed")
        network = infomapWrapper.network

        print("Building Infomap network from a NetworkX graph...")
        for e in G.edges():
            network.addLink(*e)

        print("Find communities with Infomap...")
        infomapWrapper.run()

        tree = infomapWrapper.iterTree()

        print("Found %d modules with codelength: %f" % (infomapWrapper.numTopModules(), infomapWrapper.codelength))

        communities = {}
        for node in infomapWrapper.iterLeafNodes():
            communities[node.physicalId] = node.moduleIndex()

        nx.set_node_attributes(G, name='community', values=communities)

        return infomapWrapper.numTopModules()

    def printCom(self, G, label_file):
        self.findCommunities(self.graph)
        communities = collections.defaultdict(lambda: list())
        res_comm = []
        for k, v in nx.get_node_attributes(self.graph, 'community').items():
            communities[v].append(k)
        communitie_sort = sorted(communities.values(), key=lambda b: -len(b))
        res_comm = []
        for i, c in enumerate(communitie_sort):
            res_comm.append(list(sorted(c)))
        q = self.cal_EQ(communitie_sort, G)
        res_listOfset = self.listOflist_to_listOfset(res_comm)
        onmi_value = Graph_overlap.onmi_value(res_listOfset, label_file)
        print("模块度:" + str(q))
        print("onmi:" + str(onmi_value))
        return q, onmi_value

    def cal_Pubmed(self, filename):
        file_read = open(filename, "r")
        real_label1 = dict()
        real_label = []
        for line in file_read:
            line = line.split('\t')
            real_label1[int(line[0])] = int(line[len(line) - 1].split("\n")[0])  # 从1开始要-1
        test_data_1 = sorted(real_label1.items(), key=lambda x: x[0])
        for i in range(len(test_data_1)):
            real_label.append(test_data_1[i][1])
        real_label = np.array(real_label, dtype=int)
        return real_label

    def get_reallabel(self, filename, minus=0):
        file_read = open(filename, "r")
        real_label1 = []
        for line in file_read:
            line = line.split('\t')
            real_label1.append(line[len(line) - 1].split("\n")[0])  # 从1开始要-1

        real_label = []
        for i in range(len(real_label1)):
            if real_label1[i] == '"n"':
                real_label.append(0)
            elif real_label1[i] == '"c"':
                real_label.append(1)
            elif real_label1[i] == '"l"':
                real_label.append(2)
            else:
                real_label = np.array(real_label1, dtype=int)

        if minus == 0:
            real_label = np.array(real_label)
        else:
            real_label = np.array(real_label) - 1
        return real_label

    def listOflist_to_listOfset(self, listOflist):
        setOflist = []
        listOflist = list(listOflist)
        for i in range(len(listOflist)):
            t = set()
            for j in range(len(listOflist[i])):
                t1 = int(listOflist[i][j])
                t.add(t1)
            setOflist.append(t)

        return setOflist

    def drawNetwork(self, G):
        pos = nx.spring_layout(G)
        communities = [v for k, v in nx.get_node_attributes(G, 'community').items()]
        numCommunities = max(communities) + 1
        cmapLight = colors.ListedColormap(['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed',
                                          numCommunities)
        cmapDark = colors.ListedColormap(['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed',
                                         numCommunities)

        nx.draw_networkx_edges(G, pos)

        nodeCollection = nx.draw_networkx_nodes(G,
                                                pos=pos,
                                                node_color=communities,
                                                cmap=cmapLight
                                                )
        darkColors = [cmapDark(v) for v in communities]
        nodeCollection.set_edgecolor(darkColors)

        # Draw node labels
        for n in G.nodes():
            plt.annotate(n,
                         xy=pos[n],
                         textcoords='offset points',
                         horizontalalignment='center',
                         verticalalignment='center',
                         xytext=[0, 0],
                         color=cmapDark(communities[n - 1])
                         )

        plt.axis('off')
        plt.savefig("image1.png")
        plt.show()

    def read_node_label(self, filename, skip_head=False):
        file_read = open(filename, 'r')
        X = []
        Y = []
        for line in file_read:
            if skip_head:
                line.readline()
            line = line.split('\t')
            X.append(int(line[0]) - 1)
            Y.append(int(line[1]) - 1)
        file_read.close()
        real = []
        for i in range(max(Y) + 1):
            real.append([])
        for node in range(len(Y)):
            real[Y[node]].append(node)
        return real

    def mutual_info(self, c_A, c_B, S):
        N_mA = len(c_A)
        N_mB = len(c_B)
        I_num = 0
        for i in range(len(c_A)):
            for j in range(len(c_B)):
                n_i = len(c_A[i])
                n_j = len(c_B[j])
                n_ij = len(set(c_A[i]) & set(c_B[j]))
                if n_ij == 0:
                    continue
                log_term = log((n_ij * S * 1.0) / (n_i * n_j))

                I_num += n_ij * log_term
        I_num *= -2

        I_den = 0
        for i in range(len(c_A)):
            n_i = len(c_A[i])
            I_den += n_i * log(n_i * 1.0 / S)

        for j in range(len(c_B)):
            n_j = len(c_B[j])
            I_den += n_j * log(n_j * 1.0 / S)

        I = I_num / I_den
        return I

    def cal_Q(self, partition):  # 计算Q
        m = len(self.graph.edges(None, False))
        a = []
        e = []
        for community in partition:
            t = 0.0
            for node in community:
                t += len([x for x in self.graph.neighbors(node)])
            a.append(t / (2 * m))
        for community in partition:
            t = 0.0
            for i in range(len(community)):
                for j in range(len(community)):
                    if (self.graph.has_edge(community[i], community[j])):
                        t += 1.0
            e.append(t / (2 * m))

        q = 0.0
        for ei, ai in zip(e, a):
            q += (ei - ai ** 2)
        return q

    def cal_EQ(self, cover, G):
        vertex_community = collections.defaultdict(lambda: set())
        for i, c in enumerate(cover):
            for v in c:
                vertex_community[v].add(i)

        m = 0.0
        for v in G.nodes():
            for n in G.neighbors(v):
                if v > n:
                    m += 1

        total = 0.0
        # 遍历社区
        for c in cover:
            for i in c:
                o_i = len(vertex_community[i])
                k_i = len(G[i])
                for j in c:
                    o_j = len(vertex_community[j])
                    k_j = len(G[j])
                    if i > j:
                        continue
                    t = 0.0
                    if j in G[i]:
                        t += 1.0 / (o_i * o_j)
                    t -= k_i * k_j / (2 * m * o_i * o_j)
                    if i == j:
                        total += t
                    else:
                        total += 2 * t

        return round(total / (2 * m), 4)


graph = get_G("data/youtube/youtube.graph.medium")
tool = main_tool_new.DATA(graph)
pp = tool.get_pp(3)
result_q = []
result_nmi = []
for i in range(3, 21):
    G = get_reconstructed_graph(pp, i)
    a = Infomap(G)
    label_file = "data/youtube/youtube.comm.medium"
    eq, onmi = a.printCom(graph, label_file)
    result_q.append(eq)
    result_nmi.append(onmi)
print("输出结果：")
for i in range(len(result_q)):
    print(str(result_q[i]))
for j in range(len(result_nmi)):
    print(str(result_nmi[j]))
print("最大：" + str(max(result_q)) + " " + str(max(result_nmi)))
