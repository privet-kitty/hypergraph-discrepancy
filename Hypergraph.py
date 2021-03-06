# -*- encoding: utf-8 -*-

import numpy as np
import sys
import itertools
import math

class Hypergraph:
    def __init__(self, n, m):
        self.n = n
        self.m = m
        # self.incidence = np.asarray(np.random.randint(0, 2, (n, m)))
        self.incidence = np.zeros((n, m), dtype=int)
        self.degree = 0

    def calc_degree(self):
        return np.amax(np.sum(self.incidence == 1, axis = 1))

    def delete(self, i):
        self.incidence = np.delete(self.incidence, i, 0)
        self.n -= 1
        self.degree = self.calc_degree()
        return self

    def calc_discrepancy(self, coloring):
        return np.amax(abs(np.dot(coloring, self.incidence)))

    def find_optimal_coloring(self):
        """Brute force search"""
        min_disc = self.n
        opt_coloring = np.ones(self.n, dtype=int)
        for c in itertools.product([-1, 1], repeat=self.n):
            disc = self.calc_discrepancy(c)
            if disc < min_disc:
                min_disc = disc
                opt_coloring = c
        return np.asarray(opt_coloring), min_disc

    def calc_mean_discrepancy(self):
        sum_disc = 0
        for c in itertools.product([-1, 1], repeat=self.n):
            sum_disc += self.calc_discrepancy(c)
        return sum_disc/(2**self.n)

    @staticmethod
    def roth(n):
        """Returns the set of all arimetic progressions in [n]"""
        m = 0
        for delta in range(1, n):
            for base in range(n-delta):
                len = (n-base-1)//delta
                m += len
        m += n # ROTH(n) has n unit sets.
        graph = Hypergraph(n, m)
        col = 0
        for delta in range(1, n):
            for base in range(n-delta):
                for length in range(2, (n-base-1)//delta+2):
                    for (row, idx) in zip(range(base, n, delta), range(length)):
                        graph.incidence[row, col] = 1
                    col += 1
        for base in range(n): # sequence of only a number
            graph.incidence[base, col] = 1
            col +=1
        graph.degree = graph.calc_degree()
        return graph

    @staticmethod
    def get_roth_degree(n, i):
        """Returns the degree of a given point i in ROTH(n)"""
        t = 0
        for delta in range (1, n):
            num_starting_point = math.ceil((i+1)/delta)
            num_ending_point = math.ceil((n-i)/delta)
            t += num_starting_point*num_ending_point - 1 # exclude the case: start = end
        return t+1 # add one-point case
    

    @staticmethod
    def roth_infinite(n):
        """Returns the set of all arimetic progressions without supremum nor
        infimum in Z, limited to {1, 2, ..., N}
        """
        half_n = n//2
        if (n % 2) == 1:
            m = half_n*(half_n+1) + n
        else:
            m = half_n*(half_n-1)//2 + (half_n+1)*half_n//2 + n
        
        # m = n*(n-1)//2 + 1
        graph = Hypergraph(n, m)
        col = 0
        for delta in range(1, n):
            for base in range(min(n-delta, delta)):
                for row in range(base, n, delta):
                    graph.incidence[row, col] = 1
                col += 1
        for base in range(n): # sequence of only a number
            graph.incidence[base, col] = 1
            col +=1
        graph.degree = graph.calc_degree()
        return graph

    @staticmethod
    def random(n, m):
        graph = Hypergraph(n, m)
        graph.incidence = np.asarray(np.random.randint(0, 2, (n, m)))
        graph.degree = graph.calc_degree()
        return graph

if __name__ == '__main__':
    n = int(sys.argv[1])
    print('Roth\'s Hypergraph with n =', n, ":")
    graph = Hypergraph.roth(n)
    print(graph.incidence)
    print('degree =', graph.degree)
    print('mean discrepancy =', graph.calc_mean_discrepancy())
    coloring, disc = graph.find_optimal_coloring()
    print('optimal coloring =', coloring)
    print('optimal discrepancy = ', disc)

