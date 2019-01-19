# -*- encoding: utf-8 -*-

"""
Bansal's algorithm for general hypergraph in 'Constructive Algorithms for
Discrepancy Minimization' (2010)
"""

class AlgorithmFailure(Exception):
    pass

import numpy as np
import picos as pic
import math
import random
from Hypergraph import Hypergraph

def calc_update(graph, alive_points, alpha, beta, eta, s, verbose = False):
    # Update coloring for alive variables
    incidence = graph.incidence[alive_points]
    (living_n, m) = np.shape(incidence)

    # Calculates each dangerosity of S_j
    dangerosity = np.zeros(m, dtype=int)
    max_k = beta.size-1
    for j in range(m):
        k = max(0, np.searchsorted(beta, abs(eta[j]))-1)
        if (k == max_k):
            raise AlgorithmFailure
        else:
            dangerosity[j] = k
    # print("dangerosity =", dangerosity)
    
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (living_n, living_n), vtype = 'symmetric')

    if verbose:
        print('incidence mat. of alive points =')
        print(incidence)
    
    vs = [incidence[:, j] for j in range(m)]
    vvs = [pic.new_param('vv'+str(j), np.outer(vs[j], vs[j])) for j in range(m)]
    iden = pic.new_param('I', np.identity(living_n, dtype=int))

    if verbose:
        # print("I =", iden)
        print('vs =', vs)
        print('vvs =', [np.outer(vs[j], vs[j]) for j in range(m)])

    sdp.add_constraint(iden | Xp > living_n * 0.5)
    sdp.add_list_of_constraints([vvs[j] | Xp < alpha[dangerosity[j]] for j in range(m)])
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(living_n)])
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('find', 0)
    if verbose:
        print(sdp)

    sdp.solve(verbose=0)
    v = np.linalg.cholesky(Xp.value)
    g = np.array([np.random.normal() for i in range(living_n)])
    gamma = np.zeros(graph.n)
    
    short_idx = 0
    for (i, flag) in enumerate(alive_points):
        if flag:
            gamma[i] = s * np.dot(g, v[:, short_idx])
            short_idx += 1

    if verbose:
        print('v =', v)
        print('gamma =', gamma)

    return gamma

K = 10

def sub_routine (graph, coloring, alive_points):
    n = graph.n
    m = graph.m

    a = np.sum(alive_points)
    half_a = a // 2
    print("a=", a)
    s = 1/(4*math.log2(m*n)**1.5)
    q = math.log2(2*m/a)
    print("q=", 1)
    d = 9*math.log2(20*K)
    c = 64*math.sqrt(d*(1+math.log(K)))
    print("c=", c)
    stopping_time = math.ceil(16 / s**2)
    print("stopping_time =", stopping_time)

    def beta_func (k):
        if k == 0:
            return 0
        else:
            return c * math.sqrt(a) * (q+1) * (2-1/k)

    eta = np.zeros(m) # total discrepancy incurred by S_j (allowing negative)
    ## As |eta[j]| <= n for all S_j, we have only to hold beta[0] to beta[max_k],
    ## where max_k is the smallest number that satisfies beta_func(max_k) >= n.
    max_k = 100
    for i in range(max_k):
        if beta_func(i) >= n:
            max_k = i
            break
    
    beta = np.array([beta_func (k) for k in range(max_k+1)])
    print("beta =", beta)

    alpha = np.array([d*a*(q+1)/((k+1)**5) for k in range(max_k+1)])
    print("alpha =", alpha)

    sup_barrier = 1 - 1/math.log2(m*n)
    inf_barrier = -1 + 1/math.log2(m*n)

    delta_coloring = np.zeros(n)

    for t in range(stopping_time):
        print("t=", t, "/", stopping_time)
        # print(alive_points)
        print(coloring)

        if np.sum(alive_points) <= half_a:
            return coloring
        gamma = calc_update(graph, alive_points, alpha, beta, eta, s)
        coloring += gamma
        delta_coloring += gamma
        for j in range(m):
            eta[j] = np.dot(delta_coloring, graph.incidence[:, j])
        for i in range(n):
            if abs(coloring[i]) > 1:
                raise AlgorithmFailure
            elif coloring[i] >= sup_barrier:
                alive_points[i] = False
                if random.random() <= (1 + coloring[i]) * 0.5:
                    coloring[i] = 1
                else:
                    coloring[i] = -1
            elif coloring[i] <= inf_barrier:
                alive_points[i] = False
                if random.random() <= (1 - coloring[i]) * 0.5:
                    coloring[i] = -1
                else:
                    coloring[i] = 1
    
    raise AlgorithmFailure
    
def main_routine (graph):
    n = graph.n
    m = graph.m

    x = np.zeros(n) # fractional
    l = math.ceil(math.log2(math.log2(m)))
    print("l =", l)
    alive_points = np.full(n, True)
    for i in range(0, l):
        x = sub_routine(graph, x, alive_points)

    x_int = np.copy(x)
    for i in range(0, n):
        if alive_points[i]:
            if random.random() <= (1-x[i])*0.5:
                x_int[i] = -1
            else:
                x_int[i] = 1
    return x, x_int

if __name__ == '__main__':
    # For test
    np.random.seed(0)
    graph = Hypergraph.random(6, 6)
    np.random.seed()

    print("n =", graph.n, ", m =", graph.m)
    print("Incidence matrix:\n", graph.incidence)
    print("degree=", graph.degree)

    coloring, fractional = main_routine(graph);
    print("coloring =", coloring)
    print("fractional =", fractional)
    
