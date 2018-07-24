"""
Application of Bansal-Garg framework in 'Algorithmic Discrepancy
Beyond Partial Coloring' (2016) to Beck-Fiala problem
"""

import numpy as np
import picos as pic
import math
from Hypergraph import Hypergraph

def calc_update(constraints_vectors, beta):
    (n, l) = np.shape(constraints_vectors)
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (n, n), vtype = 'symmetric')
    wws= [pic.new_param('ww'+str(i), constraints_vectors[:,i][:, np.newaxis]*constraints_vectors[:,i]) for i in range(l)]
    diagX = pic.tools.diag(pic.tools.diag_vect(Xp)/beta)
    # for i in range(dim):
    #     sdp.add_constraint(wws[i] | Xp == 0)
    sdp.add_list_of_constraints([wws[i] | Xp == 0 for i in range(l)])
    sdp.add_constraint(Xp << diagX)
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(n)], 'i')
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('max', pic.trace(Xp))
    sdp.solve(verbose=0)
    # print(Xp.value)

    return (np.linalg.cholesky(Xp.value))

def get_constraints_vector(graph, alive_points):
    alive_incidence = graph.incidence[alive_points]
    return alive_incidence[:, np.sum(alive_incidence == 1, axis = 0) >= 4*graph.degree]


def solve_beck_fiala(graph, print_per_time = 1):
    delta = 0.25
    epsilon = 0.75
    
    n = graph.n
    m = graph.m
    threshold = 1e0/n
    sup_alive = 1 - threshold
    inf_alive = -1 + threshold
    print("range of non-fixed coloring: [", sup_alive, ",", inf_alive, "]")
    gamma = 1/(n**10 * m**4 * math.log(m*n))
    # gamma=0.1
    stop_time = math.ceil((8/(epsilon*gamma*gamma))*n*math.log(n))

    coloring = np.zeros(n)
    frozen_coloring = np.zeros(n)
    alive_points = np.full(n, True)
    num_alive_points = n

    beta = (1-delta)/2
    
    print("gamma = ", gamma)
    print("T = ", stop_time)
    print("beta =", beta)

    time_with_mod = 1
    for time in range(1, 1+stop_time):
        if num_alive_points == 0:
            break;
        ws = get_constraints_vector(graph, alive_points)
        u = calc_update(ws, beta)
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
            print("t =", time, ", alive points =", num_alive_points, ",", coloring)
            time_with_mod = 0
        time_with_mod += 1

    final_coloring = np.round(coloring)
    for idx in range(n):
        if final_coloring[idx] == 0:
            final_coloring[idx] = np.random.randint(0, 2)*2-1
    print("Final coloring (float): ", coloring)
    print("Final coloring (int): ", final_coloring)

    

if __name__ == '__main__':
    graph = Hypergraph.random(2, 4)

    print("n =", graph.n, ", m =", graph.m)
    print("Incidence matrix:\n", graph.incidence)
    print("degree=", graph.degree)

    solve_beck_fiala(graph, print_per_time = 100)
