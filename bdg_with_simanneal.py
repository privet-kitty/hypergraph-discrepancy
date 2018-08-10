# -*- encoding: utf-8 -*-

"""
Bansal's algorithm with Simulated Annealing for Beck-Fiala problem.
Usage:
$ python bdg_with_simanneal.py [n]
find optimal coloring for the hypergraph ROTH(n)
"""

import numpy as np
import sys
import picos as pic
import math
from Hypergraph import Hypergraph
from simanneal import Annealer
from bansal_dadusch_garg import calc_update

def round_coloring(coloring):
    rounded_coloring = np.round(coloring).astype(int)
    rounded_coloring[rounded_coloring == 0] = 1
    return rounded_coloring

class BeckFialaProblem(Annealer):

    """
    Annealer for Beck-Fiala Problem
    """

    def __init__(self, graph):
        self.graph = graph
        n = graph.n
        m = graph.m
        self.threshold = 1e0/n
        self.sup_alive = 1 - self.threshold
        self.inf_alive = -1 + self.threshold
        print("range of non-fixed coloring: [", self.sup_alive, ",", self.inf_alive, "]")
        self.gamma = 1/(n*n*math.log(n))
        self.Tmax = -self.gamma/math.log(0.98)
        self.Tmin = -self.gamma/math.log(0.001)
        self.stop_time = math.ceil(12/(self.gamma**2))
        self.steps = self.stop_time
        self.state = np.zeros(n)
        self.alive_points = np.full(n, True) # flags
        self.num_alive_points = n
        
        print("gamma = ", self.gamma)
        print("T = ", self.stop_time)
        
        super(BeckFialaProblem, self).__init__(self.state)

    def move(self):
        if(self.num_alive_points == 0):
            self.user_exit = True
            return
        # print(self.state)
        u = calc_update(self.graph, self.alive_points, self.state)
        short_idx = 0
        r = np.random.randint(0, 2, self.num_alive_points) * 2 - 1
        for (idx, flag) in enumerate(self.alive_points):
            if flag:
                self.state[idx] += self.gamma*np.dot(r, u[short_idx])
                short_idx += 1
        for (idx, flag) in enumerate(self.alive_points):
            if flag:
                if (not (self.inf_alive <= self.state[idx]
                        and self.state[idx] <= self.sup_alive)):
                    self.alive_points[idx] = False
                    self.num_alive_points -= 1

    def energy(self):
        """Calculates the discrepancy of the coloring."""
        return self.graph.calc_discrepancy(round_coloring(self.state))

    def get_coloring(self):
        """returns (integer) coloring"""
        coloring = np.round(self.state).astype(int)
        for idx in range(self.graph.n):
            if coloring[idx] == 0:
                coloring[idx] = np.random.randint(0, 2)*2-1
        return coloring


def round_coloring(coloring):
    rounded_coloring = np.round(coloring).astype(int)
    for idx in range(coloring.size):
        if rounded_coloring[idx] == 0:
            rounded_coloring[idx] = np.random.randint(0, 2)*2-1
    return rounded_coloring

if __name__ == '__main__':

    n = int(sys.argv[1])
    graph = Hypergraph.roth(n)

    bfp = BeckFialaProblem(graph)
    bfp.copy_strategy = "slice"
    print("Tmax:", bfp.Tmax)
    print("Tmin:", bfp.Tmin)

    # start_time = time.time()
    bfp.anneal()
    state = bfp.best_state
    disc = bfp.best_energy
    # elapsed_time = time.time() - start_time

    print("n =", graph.n, ", m =", graph.m)
    print("Incidence matrix:\n", graph.incidence)
    print("degree=", graph.degree)

    print()
    print("final coloring (float):", state)
    coloring = round_coloring(state)
    print("final coloring (int):", coloring)
    print("discrepancy:", graph.calc_discrepancy(coloring))
    # print('elapsed time:', elapsed_time)
    
    opt_coloring, disc = graph.find_optimal_coloring()
    print('optimal coloring:', opt_coloring)
    print('optimal discrepancy: ', disc)

