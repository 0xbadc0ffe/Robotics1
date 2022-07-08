from sympy import *
import numpy as np
import os

os.system("cls")


def Ai(i, h):
    Ai = [
        [2*(h[i] + h[i+1]), h[i], 0],
        [h[i+2], 2*(h[i+1]+h[i+2]), h[i+1]]
    ]
    return np.asmatrix(Ai)



q_vec = [0, 2, 4, 2, 1]

v1 = 10
vn = 10

h_vec = [1,2,3,4]

N = len(q_vec)
m = len(h_vec)
assert(m == N-1)

# A1 = np.tri(m,m,1)-np.tri(m,m,0)
# A2 = np.tri(m,m,0)-np.tri(m,m,-1)
# A3 = np.tri(m,m,-1)-np.tri(m,m,-2)
A = np.zeros((m-1,m-1))
# print(A1)
# print(A2)
# print(A3)

A[0,0] = 2*(h_vec[0] + h_vec[1])

for i in range(m-2):
    A[i,i+1]=h_vec[i]
    A[i+1,i+1]=2*(h_vec[i+1] + h_vec[i+2])
    A[i+1,i]=h_vec[i+2]

print("\nA:")
print(A)



b = np.zeros(m-1)
for i in range(m-2):
    b[i] = 3/(h_vec[i]*h_vec[i+1])*(h_vec[i]**2*(q_vec[i+2]-q_vec[i+1]) + h_vec[i+1]**2*(q_vec[i+1]-q_vec[i]))
b[0] -= h_vec[1]*v1
b[-1] -= h_vec[-1]*vn

print("\nb:")
print(b)



v = []
for k in range(N-2):
    v.append(symbols(f'v{k+2}'))
#v2, v3, v4 = symbols('v2, v3, v4')   # must be N-2
system = Matrix(A), Matrix(b)
s = list(linsolve(system, *v))[0]
v = list(s)
print("\nvk:")
v = [v1] + v + [vn]
print(v)

a = []
for k in range(N-2):
    a.append(symbols(f'a{k+1}2'))
    a.append(symbols(f'a{k+1}3'))

for k in range(m):
    A = [
        [h_vec[i]**2, h_vec[i]**3],
        [2*h_vec[i], 3*h_vec[i]**2]
    ]
    b = [q_vec[k+1]-q_vec[k]-v[k]*h_vec[k], v[k+1]-v[k]]
    system = Matrix(A), Matrix(b)
    s = list(linsolve(system, *a[2*k:2*k+2]))[0]
    print(f"a{k+1}2, a{k+1}3:  ",s)
    #a[f"a{k}"]