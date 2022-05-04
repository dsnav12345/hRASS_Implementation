from docplex.mp.model import Model
from graphviz import Digraph


def vGDA(C, Ct, m, e):
    r = len(C)
    wx = [1]
    d = 0
    T = 0
    while True:
        d += 1
        wx.insert(0, wx[0]*m)
        T = int(Ct*wx[0])
        if wx[0]*e*Ct > abs(wx[0]*Ct - T):
            break
        T += 1
        if wx[0]*e*Ct > abs(wx[0]*Ct - T):
            break

    B = []
    for c in C:
        B.append(int(c*wx[0]))

    print(B, T, d)
    W = []
    for i in range(r):
        W.append([])
        W[i].append(B[i])
        for j in range(d):
            W[i].append(W[i][j]/m)

    obj = Model('vGDA')

    X = []

    for i in range(r):
        X.append([])
        for j in range(d+1):
            X[i].append(obj.integer_var(lb=0, ub=m))

    obj.add_constraint(sum(X[i][j]*W[i][j] for i in range(r) for j in range(d+1))<=T+2, ctname='c1')
    obj.add_constraint(sum(X[i][j]*W[i][j] for i in range(r) for j in range(d+1))>=T-2, ctname='c2')
    obj.add_constraint(sum(X[i][j]*wx[j] for i in range(r) for j in range(d+1))<=wx[0], ctname='c3')

    cpx = obj.get_cplex()
    cpx.populate_solution_pool()
    n = cpx.solution.pool.get_num()

    ans = []
    a_m = (d+1)*r

    for idx in range(n):
        ls = cpx.solution.pool.get_values(idx)
        t_sol = []
        for i in range(r):
            t_sol.append([])
            for j in range(d+1):
                t_sol[i].append(int(ls[i*(d+1)+j]))

        n_m = 0
        prev = 0
        for i in range(d, 0, -1):
            prev += sum(t_sol[j][i] for j in range(r-1))
            prev = prev//m + int(bool(prev%m))
            n_m += prev
        if a_m > n_m:
            ans = t_sol
            a_m = n_m

    return [a_m, ans]

def build_graph(m, X):
    g = Digraph()
    r = len(X)
    d = len(X[0])
    prev = []
    idx = 0
    nm = 1
    while sum(X[i][d-1] for i in range(r)) == 0:
        d -= 1

    for i in range(d-1, 0, -1):
        tem = m
        g.node(str(idx), 'mix'+str(nm))
        nm += 1
        curr = idx
        for j in prev:
            g.edge(str(j), str(idx), label='1')
            tem -= 1
            if tem == 0:
                tem = m
                g.node(str(idx), 'mix'+str(nm))
                nm += 1
                curr = idx
                idx += 1
        prev.clear()
        prev.append(idx)
        idx += 1
        for j in range(r):
            if tem == 0:
                tem = m
                g.node(str(idx))
                curr = idx
                idx += 1
            if X[j][i] > 0:
                g.node(str(idx), 'A{}'.format(j+1))
                g.edge(str(idx), str(curr), label=str(min(X[j][i],tem)))
                idx += 1
                tem -= X[j][i]
                if tem < 0:
                    g.node(str(idx), 'mix'+str(nm))
                    nm += 1
                    curr = idx
                    prev.append(idx)
                    idx += 1
                    g.node(str(idx), 'A{}'.format(j+1))
                    g.edge(str(idx), str(curr), str(0-tem))
                    idx += 1
                    tem += m
        if 0 < tem < m:
            g.node(str(idx), 'B')
            g.edge(str(idx), str(curr), label=str(tem))
            idx += 1
    g.view()


Ct = float(input('Final concentration : '))
C = list(map(float, input('Reagent Concentrations : ').split()))
M = int(input('Mixer size : '))
e = float(input('Error tolerance : '))

f_ans = None
o_nm = len(C)*M
o_m = 0
for m in range(M, 1, -1):
    tem = vGDA(C, Ct, m, e)
    if tem[0] < o_nm:
        o_nm = tem[0]
        f_ans = tem[1]
        o_m = m

print('mixer size = {}\nxij={}'.format(o_m, f_ans))
build_graph(o_m, f_ans)
