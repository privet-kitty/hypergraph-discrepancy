
import numpy as np
import picos as pic
import math


def get_matrix_one(m, n, i, j):
    mat = np.zeros((m, n))
    mat[i, j] = 1
    return mat

def calc_update():

    delta = 0.25
    beta = (1-delta)/2
    
    dim = 4
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (dim, dim), vtype = 'symmetric')
    ws = np.array([[0, 1, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [1, 0, 0, 1]])
    ws = np.transpose(ws)
    wws= [pic.new_param('ww'+str(i), ws[:,i][:, np.newaxis]*ws[:,i]) for i in range(dim)]
    print(wws[0])
    # ones = [get_matrix_one(dim, dim, i, i) for i in range(dim)]
    diagX = pic.tools.diag(pic.tools.diag_vect(Xp)/beta)
    for i in range(dim):
        sdp.add_constraint(wws[i] | Xp == 0)
    # sdp.add_list_of_constraints([wws[i] | Xp == 0 for i in range(dim)])
    sdp.add_constraint(Xp << diagX)
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(dim)], 'i')
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('max', pic.trace(Xp))
    sdp.solve()
    print(sdp)
    print(Xp)

    return(Xp)


def calc_update2():

    delta = 0.25
    beta = (1-delta)/2
    
    dim = 4
    sdp = pic.Problem()
    Xp = sdp.add_variable("X", (dim, dim), vtype = 'symmetric')
    ws = np.array([[0, 1, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [1, 0, 0, 1]])
    wws= [pic.new_param('ww'+str(i), ws[i][:, np.newaxis]*ws[i]) for i in range(dim)]
    print(wws[0])
    # ones = [get_matrix_one(dim, dim, i, i) for i in range(dim)]
    diagX = pic.tools.diag(pic.tools.diag_vect(Xp)/beta)
    for i in range(dim):
        sdp.add_constraint(wws[i] | Xp == 0)
    # sdp.add_list_of_constraints([wws[i] | Xp == 0 for i in range(dim)])
    sdp.add_constraint(Xp << diagX)
    sdp.add_list_of_constraints([Xp[i, i] < 1 for i in range(dim)], 'i')
    sdp.add_constraint(Xp >> 0)
    sdp.set_objective('max', pic.trace(Xp))
    sdp.solve()
    print(sdp)
    print(Xp)

    return(Xp)

calc_update()
calc_update2()
