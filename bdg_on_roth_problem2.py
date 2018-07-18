# -*- encoding: utf-8 -*-

#
# Efficiency test of Bansal's algorithm for Roth's problem
#
# usage: python bdg_on_roth_problem2.py [end_n] [sample]
# => investigate ROTH[2] to ROTH[end_n-1] with sample = [sample]

import sys
import time as tm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pickle
from Hypergraph import Hypergraph
from bansal_dadusch_garg import *

def solve_roth(n):
    graph = Hypergraph.roth(n)
    print('n =', n, ':')
    print(graph.incidence)
    start = tm.time()
    coloring, num_steps = solve(graph, print_per_time = 100)
    elapsed_time = tm.time() - start
    return coloring, num_steps, elapsed_time

def solve_roth_bf(n):
    # brute force
    graph = Hypergraph.roth(n)
    print('n =', n, ':')
    print(graph.incidence)
    start = tm.time()
    coloring = graph.find_optimal_coloring()
    elapsed_time = tm.time() - start
    return coloring, elapsed_time

def measure_time(n, sample = 10):
    result = np.asarray([solve_roth(n) for i in range(sample)])
    num_steps_median = np.median(result[:, 1])
    elapsed_time_median = np.median(result[:, 2])
    return num_steps_median, elapsed_time_median

def measure_time_bf(n, sample = 10):
    result = np.asarray([solve_roth_bf(n) for i in range(sample)])
    elapsed_time_median = np.median(result[:, 1])
    return elapsed_time_median


sample = int(sys.argv[2])
nrange = np.arange(2, int(sys.argv[1]))
result = np.asarray([measure_time(n, sample = sample) for n in nrange])
elapsed_time_bf = np.asarray([measure_time_bf(n, sample = sample) for n in nrange])
num_steps = result[:, 0]
elapsed_time = result[:, 1]
print(elapsed_time)
print(elapsed_time_bf)
print(num_steps)

with open("elapsed_time.pickle", "wb") as f:
    pickle.dump(elapsed_time, f)
with open("elapsed_time_bf.pickle", "wb") as f:
    pickle.dump(elapsed_time_bf, f)
with open("num_steps_lst.pickle", "wb") as f:
    pickle.dump(num_steps, f)

plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))

plt.plot(nrange, elapsed_time, label="Bansal")
plt.plot(nrange, elapsed_time_bf, label="brute force")
plt.title('Elapsed time')
plt.xlabel('n')
plt.ylabel('sec.')
plt.legend()
plt.savefig('elapsed_time.png')
plt.show()

plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
plt.plot(nrange, num_steps)
plt.title('Number of solve(SDP)')
plt.xlabel('n')
plt.ylabel('number')
plt.savefig('num_steps.png')
plt.show()
