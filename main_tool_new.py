# coding=utf-8
import networkx as nx
import numpy as np
# import cupy as np
from scipy.sparse import csr_matrix
# from cupyx.scipy.sparse import csr_matrix
import Graph
import igraph as ig
# from matplotlib import pyplot as plt
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "2"


def matrix_point_multiply(a_matrix, a_vector):
    return np.einsum('ij,i->ij', a_matrix.T, a_vector).T


class GRAPH(object):
    G_nx = []
    G_ig = []
    G = []
    G1 = nx.Graph()
    G1.add_edge(1, 2)
    G1.add_edge(1, 3)
    G1.add_edge(3, 2)
    G1.add_edge(4, 3)
    G1.add_edge(4, 5)
    G1.add_edge(4, 6)
    G1.add_edge(5, 7)

    def __init__(self, network_filename):
        self.G_nx, self.G_ig = Graph.my_graph(network_filename)
        self.G = self.G_nx

    @staticmethod
    def nx2ig(nx_graph):
        d = nx.to_pandas_edgelist(nx_graph).values
        d = d - 1
        newG = ig.Graph(d)
        return newG


class DATA(object):
    N = 0
    AP = []
    FAP = []
    ap_current_number = 1
    A = csr_matrix([[]])
    B = csr_matrix([[]])
    BT = B.T
    BP = []
    PP = [[]]  #
    pp_current_number = 1

    def __init__(self, g):

        self.N = len(g)

        self.A = csr_matrix(np.array(nx.adjacency_matrix(g).todense()))
        en = np.ones(self.N)
        du = self.A.dot(en)
        self.Du = du
        np.where(du == 0, 1, du)
        self.BT = csr_matrix(self.A.T / du)
        self.B = self.BT.T

        self.AP = [[] for _ in range(500)]
        self.AP[1] = self.B.toarray()

        self.BP = [[] for _ in range(500)]
        self.__compute_bp()

        self.FAP = [[] for _ in range(500)]
        self.FAP[1] = self.AP[1]

        self.PP = [[] for _ in range(500)]
        self.PP[1] = self.AP[1]

    def __compute_bp(self):
        self.BP[self.ap_current_number] = self.AP[self.ap_current_number].diagonal()

    def __compute_ap(self, n):
        the_last = self.ap_current_number
        for _ in range(n - the_last):
            self.AP[self.ap_current_number + 1] = self.BT.dot(self.AP[self.ap_current_number].T).T
            self.ap_current_number += 1
            self.__compute_bp()

    def __compute_fap(self, n):
        the_last = self.pp_current_number
        for j in range(the_last + 1, n + 1):
            ap = self.get_ap(j)
            temp = np.array([[0 for _ in range(self.N)] for _ in range(self.N)])
            for k in range(1, j - 1):
                temp = matrix_point_multiply(self.FAP[k], self.BP[j - k]) + temp
            self.FAP[j] = ap - temp
        del temp

    def __compute_pp(self, n):
        the_last = self.pp_current_number
        for i in range(the_last + 1, n + 1):  # type: int
            last_pp1 = self.PP[i - 1]
            fap = self.__get_fap(i)
            self.pp_current_number += 1
            self.PP[i] = last_pp1 + fap

    def get_a(self):
        return self.A.toarray()

    def get_b(self):
        return self.B.toarray()

    def get_ap(self, n):
        if n > self.ap_current_number:
            self.__compute_ap(n)
        return self.AP[n]

    def get_bp(self, n):
        if n > self.ap_current_number:
            self.__compute_ap(n)
        return self.BP[n]

    def __get_fap(self, n):
        if n > self.pp_current_number:
            self.__compute_fap(n)
        return self.FAP[n]

    def get_pp(self, n):
        if n > self.pp_current_number:
            self.__compute_pp(n)
        return self.PP[n]

    def get_mp(self, n):
        if n > self.ap_current_number:
            self.__compute_ap(n)
        return self.AP[n].T * self.AP[n]

    def get_tp(self, n):
        if n > self.pp_current_number:
            self.__compute_pp(n)
        return self.PP[n].T * self.PP[n]


class GRModel:
    def __init__(self, nx_graph):
        nx_G = nx_graph
        self.N = len(nx_G)

    def get_reconstructed_graph(self, resourceMatrix, num_edge=1):
        G1 = nx.Graph()
        for i in range(self.N):
            max_key = [0] * num_edge
            max_value = [0] * num_edge
            for j in range(num_edge):
                max_key[j] = np.argmax(resourceMatrix[i])
                max_value[j] = resourceMatrix[i][max_key[j]]
                resourceMatrix[i][max_key[j]] = 0
            for j in range(num_edge):
                resourceMatrix[i][max_key[j]] = max_value[j]
                G1.add_edge(i + 1, max_key[j] + 1)
        return G1
