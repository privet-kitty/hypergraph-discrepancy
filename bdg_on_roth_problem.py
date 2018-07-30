# -*- encoding: utf-8 -*-

"""
Basic test of Bansal's algorithm for Roth's problem
Usage: 
$ python bdg_on_roth_problem.py [n]
"""

import sys
from Hypergraph import Hypergraph
from bansal_dadusch_garg import *

n = int(sys.argv[1])
graph = Hypergraph.roth(n)

print("n =", graph.n, ", m =", graph.m)
print("Incidence matrix:\n", graph.incidence)
print("degree=", graph.degree)

coloring, time = solve(graph, print_per_time = 100)
print("final discrepancy =", graph.calc_discrepancy(coloring))

opt_coloring, disc = graph.find_optimal_coloring()
print('optimal coloring =', opt_coloring)
print('optimal discrepancy = ', disc)

