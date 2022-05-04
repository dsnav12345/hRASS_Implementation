from docplex.mp.model import Model
from graphviz import Digraph


def vGDA(B, T, m, d):
    r = len(B)
    wx = [1]
    for i in range(d):
        wx.insert(0, wx[0]*m)

    # print(B, T, d)
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
            X[i].append(obj.integer_var(lb=0, ub=m-1))

    obj.add_constraint(sum(X[i][j]*W[i][j] for i in range(r) for j in range(d+1))==T, ctname='c1')
    obj.add_constraint(sum(X[i][j]*wx[j] for i in range(r) for j in range(d+1))<=wx[0], ctname='c3')
    # obj.set_objective('min', (sum(X[i][j]*W[i][j] for i in range(r) for j in range(d+1))-T)*(sum(X[i][j]*W[i][j] for i in range(r) for j in range(d+1))-T))

    cpx = obj.get_cplex()
    cpx.populate_solution_pool()
    n = cpx.solution.pool.get_num()

    ans = []
    a_m = (d+1)*r

    print(n, end=' ')
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
            prev += sum(t_sol[j][i] for j in range(r))
            prev = prev//m + int(bool(prev%m))
            n_m += prev
        if a_m > n_m:
            ans = t_sol
            a_m = n_m

    return [a_m, ans]


def hRASS(B, T, M, d):
    f_ans = None
    o_nm = len(B)*M
    o_m = 0
    for m in range(M, 1, -1):
        tem = vGDA(B, T, m, d)
        if tem[0] < o_nm:
            o_nm = tem[0]
            f_ans = tem[1]
            o_m = m
    return [o_m, f_ans]


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
                g.node(str(idx), 'A{}'.format(j))
                g.edge(str(idx), str(curr), label=str(min(X[j][i],tem)))
                idx += 1
                tem -= X[j][i]
                if tem < 0:
                    g.node(str(idx), 'mix'+str(nm))
                    nm += 1
                    curr = idx
                    prev.append(idx)
                    idx += 1
                    g.node(str(idx), 'A{}'.format(j))
                    g.edge(str(idx), str(curr), str(0-tem))
                    idx += 1
                    tem += m
        if 0 < tem < m:
            g.node(str(idx), 'B')
            g.edge(str(idx), str(curr), label=str(tem))
            idx += 1
    g.view()


f = open('input.txt')

ls = list(map(int, f.readline().split()))
M = ls[0]
d = ls[1]
R = ls[2]
B = []
for i in range(3, R+3):
    B.append(ls[i])

T = ls[R+3]

ans = vGDA(B, T, M, d)
print(ans[0], ans[1])
build_graph(M, ans[1])

# f = open('../Tests/6_3_3.txt')
#
# while True:
#     ls = list(map(int, f.readline().split()))
#     if len(ls)==0:
#         break
#     M = ls[0]
#     d = ls[1]
#     R = ls[2]
#     B = []
#     for i in range(3, R+3):
#         B.append(ls[i])
#
#     T = ls[R+3]
#
#     ans = vGDA(B, T, M, d)
#     print(ans[0], ans[1])
