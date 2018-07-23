import math
import sys
import numpy as np
from Hypergraph import Hypergraph
from simanneal import Annealer

"""
Let ROTH(n) be the set of all arithmetic progressions in {1, 2, ...,
n}.  It finds a coloring of small discprepancy for ROTH(n) by
simulated annealing.

Usage
-----
$ python3 vbp_simmanneal.py [n]
"""

class BeckFialaProblem(Annealer):

    """
    Annealer for Beck-Fiala Problem
    """

    def __init__(self, state, graph):
        self.graph = graph
        super(BeckFialaProblem, self).__init__(state)

    def move(self):
        """Inverses a color"""
        idx = np.random.randint(0, len(self.state)-1)
        self.state[idx] = -self.state[idx]

    def energy(self):
        """Calculates the discrepancy of the coloring."""
        return self.graph.calc_discrepancy(self.state)



if __name__ == '__main__':

    n = int(sys.argv[1])
    graph = Hypergraph.roth(n)
    init_coloring = np.asarray([np.random.randint(0, 2)*2-1 for i in range(n)])
    print("init. coloring:", init_coloring) 

    bfp = BeckFialaProblem(init_coloring, graph)

    bfp.steps = 100000
    bfp.copy_strategy = "slice"
    state, disc = bfp.anneal()

    print("n =", graph.n, ", m =", graph.m)
    print("Incidence matrix:\n", graph.incidence)
    print("degree=", graph.degree)

    print()
    print("coloring by simanneal:", state)
    print("discrepancy by simanneal:", disc)
    
    opt_coloring, disc = graph.find_optimal_coloring()
    print('optimal coloring =', opt_coloring)
    print('optimal discrepancy = ', disc)


