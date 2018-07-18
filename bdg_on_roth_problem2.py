# -*- encoding: utf-8 -*-

#
# Efficiency test of Bansal's algorithm for Roth's problem
#

import sys
import time as tm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

def measure_time(n, sample = 10):
    result = np.asarray([solve_roth(n) for i in range(sample)])
    elapsed_time_median = np.median(result[:, 2])
    num_steps_median = np.median(result[:, 1])
    return num_steps_median, elapsed_time_median


nrange = np.arange(2, 5)
result = np.asarray([measure_time(n, sample = 2) for n in nrange])
num_steps_lst = result[:, 0]
elapsed_time_lst = result[:, 1]
print(elapsed_time_lst)
print(num_steps_lst)

with open("elapsed_time.dat", "w") as f:
    for (n, time) in zip(nrange, elapsed_time_lst):
        print(n, time, file=f)

with open("num_steps_lst.dat", "w") as f:
    for (n, steps) in zip(nrange, num_steps_lst):
        print(n, steps, file=f)

plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))

plt.plot(nrange, elapsed_time_lst)
plt.title('Elapsed time')
plt.xlabel('n')
plt.ylabel('sec.')
plt.savefig('elapsed_time.png')
plt.show()

plt.gca().get_xaxis().set_major_locator(ticker.MaxNLocator(integer=True))
plt.plot(nrange, num_steps_lst)
plt.title('Number of solving SDPs')
plt.xlabel('n')
plt.ylabel('Number')
plt.savefig('num_steps.png')
plt.show()
