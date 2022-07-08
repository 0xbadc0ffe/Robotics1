from sympy import *
from sympy.solvers.solveset import linsolve


# RESOLVING THIS SYS OF EQ.
#
# x + y + z = 1 
# x + y + 2z = 3
#


x, y, z = symbols('x, y, z')

#List of Equations Form:

s = linsolve([x + y + z - 1, x + y + 2*z - 3 ], (x, y, z))
print(s)

#Augmented Matrix Form:
s = linsolve(Matrix(([1, 1, 1, 1], [1, 1, 2, 3])), (x, y, z))
print(s)

#A*x = b Form
M = Matrix(((1, 1, 1, 1), (1, 1, 2, 3)))
system = A, b = M[:, :-1], M[:, -1]
s = linsolve(system, x, y, z)
print(s)


# TEST

a1, b1, a2, b2 = symbols('a1, b1, a2, b2')
th_v, th_1, th_2, vT1, vT2  = symbols('θv, θ1, θ2, vT1, vT2')
s = linsolve([a1+b1-th_v+th_1, 3*a1+2*b1-vT1, -a2+b2-th_v+th_2, 3*a2-2*b2-vT2], (a1, b1, a2, b2))
print(s)