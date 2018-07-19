# -*- encoding: utf-8 -*-

#
# Accuracy test of Bansal's algorithm for Roth's problem
#
# usage: python bdg_on_roth_problem3.py [n] [sample]

import sys
import time as tm
import numpy as np
from Hypergraph import Hypergraph
import bansal_dadusch_garg as bdg
from scipy import stats

n = int(sys.argv[1])
sample = int(sys.argv[2])
graph = Hypergraph.roth(n)

def solve_and_print(i, graph):
    print('[Trial ', i, ']', sep='')
    return bdg.solve(graph, print_per_time = 100)[0]
    
colorings = np.asarray([solve_and_print(i, graph) for i in range(sample)])
discs = [graph.calc_discrepancy(colorings[i]) for i in range(sample)]

mean = np.mean(discs)
mean_random = graph.calc_mean_discrepancy()
sd = np.std(discs)
t = (mean-mean_random)/(sd/np.sqrt(sample))
df = sample-1
print('')
print('[Test for ROTH(', n, ')]', sep='')
print('sample =', sample)
print('discrepancies by Bansal =', discs)
print('mean disc. by Bansal =', mean)
print('mean disc. by random coloring =', graph.calc_mean_discrepancy())
print('H0 (null hypothesis): mean disc. by Bansal >= mean disc. of random coloring')
print('t stat. =', t)
ppf = stats.t.ppf(0.01, df)
print('ppf of t-distribution with p<0.01, 1-tail, df=', df, ': ', ppf, sep='')
if t < ppf:
    print('H0 is rejected.')
else:
    print('H0 is not rejected.')
