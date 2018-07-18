# -*- encoding: utf-8 -*-
#
# Checks if PICOS returns the same result for an SDP and its standard
# form.
#

import numpy as np
import picos as pic
import math

class Hypergraph:
    def __init__(self, n, m):
        self.n = n
        self.m = m
        self.incidence = np.asarray(np.random.randint(0, 2, (n, m)))

    def degree(self):
        return np.amax(np.sum(self.incidence == 1, axis = 1))

    def delete(self, i):
        self.incidence = np.delete(self.incidence, i, 0)
        self.n -= 1
        return self

def calc_update(graph):
    incidence = graph.incidence
    print(incidence)
    (n, l) = np.shape(incidence)
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (n, n), vtype = 'symmetric')
    a = 6
    degree = graph.degree()
    big_set_flags = (np.sum(incidence == 1, axis = 0) >= a*degree)
    coloring = np.asarray(np.random.rand(n))*0.4-0.2

    vs = [incidence[:, i] for i in range(l)]
    xs = [vs[i]*coloring for i in range(l)]
    vvs = [pic.new_param('vv'+str(i), np.outer(vs[i], vs[i])) for i in range(l)]
    iden = pic.new_param('I', np.identity(n))
    xxs = [pic.new_param('xx'+str(i), np.outer(xs[i], xs[i])) for i in range(l)]

    # print('vs =')
    # print(vs)
    # print('xs =')
    # print(xs)
    # print('vvs =')
    # print([np.outer(vs[i], vs[i]) -2*np.identity(n) for i in range(l)])
    # print('xxs =')
    # print([np.outer(xs[i], xs[i]) -2*np.identity(n) for i in range(l)])

    sdp.add_list_of_constraints([vvs[i] | Xp == 0 for i in range(l) if big_set_flags[i]])
    sdp.add_list_of_constraints([(vvs[i] - 2*iden) | Xp < 0 for i in range(l) if not big_set_flags[i]])
    sdp.add_list_of_constraints([(xxs[i] - 2*iden) | Xp < 0 for i in range(l) if not big_set_flags[i]])
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(n)])
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('max', pic.trace(Xp))
    sdp.solve()
    print(Xp.value)
    U = np.linalg.cholesky(Xp.value)
    print(U)
    
    return U, np.asarray(Xp.value)

def make_one(row, col, coord):
    mat = np.zeros((row, col), dtype=int)
    mat[coord][coord] = 1
    return mat

def calc_update_standard(graph):
    print('[standard form]')
    incidence = graph.incidence
    print(incidence)
    (n, l) = np.shape(incidence)
    a = 6
    degree = graph.degree()
    big_set_flags = (np.sum(incidence == 1, axis = 0) >= a*degree)
    num_big = np.sum(big_set_flags)
    num_small = l - num_big
    num_pad = num_big + 2*num_small + n
    total_size = num_pad + n
    padding = ((0, num_pad), (0, num_pad))
    coloring = np.asarray(np.random.rand(n))*0.4-0.2

    print('num_big =')
    print(num_big)
    print('num_pad =')
    print(num_pad)
    print('total_size =')
    print(total_size)
          
    sdp = pic.Problem()
    Zp = sdp.add_variable("X", (total_size, total_size), vtype = 'symmetric')
    
    vs = [incidence[:, i] for i in range(l)]
    xs = [vs[i]*coloring for i in range(l)]
    vvs = [np.pad(np.outer(vs[i], vs[i]), padding, mode='constant') for i in range(l)]
    xxs = [np.pad(np.outer(xs[i], xs[i]), padding, mode='constant') for i in range(l)]
    extra= [make_one(total_size, total_size, n+i) for i in range(num_pad)]
    single_one = [make_one(total_size, total_size, i) for i in range(n)]
    
    vvs_p = [pic.new_param('vv'+str(i), vvs[i]) for i in range(l)]
    iden_p = pic.new_param('I', np.pad(np.identity(n), padding, mode='constant'))
    xxs_p = [pic.new_param('xx'+str(i), xxs[i]) for i in range(l)]
    extra_p = [pic.new_param('extra-'+str(i), extra[i]) for i in range(num_pad)]
    single_one_p = [pic.new_param('single-one-'+str(i), single_one[i]) for i in range(n)] 

    # print('vs =')
    # print(vs)
    # print('xs =')
    # print(xs)
    # print('vvs =')
    # print(vvs)
    # print('single_one =')
    # print(single_one)
    # print('extra_part =')
    # print(extra)

    sdp.add_list_of_constraints([vvs_p[i] | Zp == 0 for i in range(l) if big_set_flags[i]])
    extra_idx=0
    for i in range(l):
        if not big_set_flags[i]:
            sdp.add_constraint((vvs_p[i] - 2*iden_p + extra_p[extra_idx]) | Zp == 0)
            extra_idx += 1
    for i in range(l):
        if not big_set_flags[i]:
            sdp.add_constraint((xxs_p[i] - 2*iden_p + extra_p[extra_idx]) | Zp == 0)
            extra_idx += 1
    for i in range(n):
        sdp.add_constraint((single_one_p[i] + extra_p[extra_idx]) | Zp == 1)
        # sdp.add_constraint(Zp[i, i] + Zp[n+extra_idx, n+extra_idx] == 1)
        extra_idx += 1
    sdp.add_constraint(Zp >> 0)
    sdp.set_objective('max', iden_p | Zp)

    print(sdp)

    sdp.solve()

    X = Zp.value[0:n, 0:n]
    print(X)
    U = np.linalg.cholesky(X)
    print(U)
    
    return U, np.asarray(X)

graph = Hypergraph(15, 20)
U1, X1 = calc_update(graph)
U2, X2 = calc_update_standard(graph)
print(U1-U2)
