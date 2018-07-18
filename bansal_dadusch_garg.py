# -*- encoding: utf-8 -*-
#
# Bansal's algorithm for Beck-Fiala problem in 'An Algorithm for
# Komlós Conjecture Matching Banaszczyk's Bound' (2016)
#

import numpy as np
import picos as pic
import math
from Hypergraph import Hypergraph

def calc_update(graph, alive_points, coloring, verbose = False):
    # Updates alive variables
    incidence = graph.incidence[alive_points]
    (remaining_n, l) = np.shape(incidence)
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (remaining_n, remaining_n), vtype = 'symmetric')
    a = 6
    big_set_flags = (np.sum(incidence == 1, axis = 0) >= a*graph.degree)
    alive_coloring = coloring[alive_points]
    if verbose:
        print('incidence mat. of alive points =')
        print(incidence)
    
    vs = [incidence[:, i] for i in range(l)]
    xs = [vs[i]*alive_coloring for i in range(l)]
    vvs = [pic.new_param('vv'+str(i), np.outer(vs[i], vs[i])) for i in range(l)]
    iden = pic.new_param('I', np.identity(remaining_n))
    xxs = [pic.new_param('xx'+str(i), np.outer(xs[i], xs[i])) for i in range(l)]

    if verbose:
        print('vs =')
        print(vs)
        print('xs =')
        print(xs)
        print('vvs =')
        print([np.outer(vs[i], vs[i]) for i in range(l)])
        print('xxs =')
        print([np.outer(xs[i], xs[i]) for i in range(l)])
    
    sdp.add_list_of_constraints([vvs[i] | Xp == 0 for i in range(l) if big_set_flags[i]])
    sdp.add_list_of_constraints([(vvs[i] - 2*iden) | Xp < 0 for i in range(l) if not big_set_flags[i]])
    sdp.add_list_of_constraints([(xxs[i] - 2*iden) | Xp < 0 for i in range(l) if not big_set_flags[i]])
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(remaining_n)])
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('max', pic.trace(Xp))

    sdp.solve(verbose=0)
    U = np.linalg.cholesky(Xp.value)

    if verbose:
        print(sdp)
        print('U =')
        print(U)

    return (np.linalg.cholesky(Xp.value))


def solve(graph, print_per_time = 1):
    n = graph.n
    m = graph.m

    threshold = 1e0/n
    sup_alive = 1 - threshold
    inf_alive = -1 + threshold
    print("range of non-fixed coloring: [", sup_alive, ",", inf_alive, "]")
    gamma = 1/(n*n*math.log(n))
    stop_time = math.ceil(12/(gamma**2))

    coloring = np.zeros(n)
    alive_points = np.full(n, True) # flags
    num_alive_points = n

    print("gamma = ", gamma)
    print("T = ", stop_time)

    time_with_mod = 1
    for time in range(1, stop_time+1):
        if num_alive_points == 0:
            break;
        u = calc_update(graph, alive_points, coloring)
        short_idx = 0
        r = np.random.randint(0, 2, num_alive_points) * 2 - 1
        for (idx, flag) in enumerate(alive_points):
            if flag:
                coloring[idx] += gamma*np.dot(r, u[short_idx])
                short_idx += 1
        for (idx, flag) in enumerate(alive_points):
            if flag:
                if (not (inf_alive <= coloring[idx] and coloring[idx] <= sup_alive)):
                    alive_points[idx] = False
                    num_alive_points -= 1
        if time_with_mod == print_per_time:
            print("t =", time, ", alive points =", num_alive_points, ", coloring =", coloring)
            time_with_mod = 0
        time_with_mod += 1

    final_coloring = np.round(coloring).astype(int)
    for idx in range(n):
        if final_coloring[idx] == 0:
            final_coloring[idx] = np.random.randint(0, 2)*2-1
    print("Stop time:", time)
    print("Final coloring (float):", coloring)
    print("Final coloring (int):", final_coloring)
    return final_coloring, time


if __name__ == '__main__':
    # For test
    graph = Hypergraph.random(6, 6)

    print("n =", graph.n, ", m =", graph.m)
    print("Incidence matrix:\n", graph.incidence)
    print("degree=", graph.degree)

    solve(graph, print_per_time = 100)
