from re import A
from sympy import *
import os, platform

if platform.system() == 'Windows':
    CLEAR_STR = "cls" 
else:
    CLEAR_STR = "clear"


### ES1
os.system(CLEAR_STR)

print("\n\n###########  ES1")


phi = Symbol("phi")
theta = Symbol("theta")
psi = Symbol("psi")

T = [
    [0, -sin(phi),  cos(phi)*cos(theta)],
    [0, cos(phi),   sin(phi)*cos(theta)],
    [1,         0,          -sin(theta)]
]

T = Matrix(T)
print("\nT:")
pprint(T)

print("\ndet(T):")
print(det(T).simplify())
print()
Ts1 = T.subs(theta,pi/2)
Ts2 = T.subs(theta,-pi/2)
pprint(Ts1)
pprint(Ts2)
print()
print("\nR{T}:")
pprint(Ts1.columnspace())
print("\nKer{T}:")
pprint(Ts1.nullspace())
print()
print("\nT^-1:")
pprint(T.inv())
#pprint(det(T.inv()).simplify())

input("\n\npress Enter")








##### ES2
os.system(CLEAR_STR)

print("\n\n###########  ES2")

#print((cos(theta+phi)*cos(theta) + sin(theta)*sin(theta+phi)).simplify())



p = [0.3,-0.3,0.7]
M=0.5
N=0.5
L=0.5
q3 = asin((p[2]-M)/N)
N2 = N*cos(q3)


#q2 = acos((p[0]**2+p[1]**2-L**2-N**2*cos(q3)**2)/(2*L*cos(q3)*N))

c2 = (p[0]**2+p[1]**2-L**2-N**2*cos(q3)**2)/(2*L*cos(q3)*N)
s2 = +sqrt(1-c2**2)  #-
q2 = atan2(s2,c2)


s1 = ((L + N2*c2)*p[1] -N2*s2*p[0])/(p[0]**2+p[1]**2)
c1 = ((L + N2*c2)*p[0] +N2*s2*p[1])/(p[0]**2+p[1]**2)
q1 = atan2(s1,c1)

# print("\n############## q1,q2,q3")
# print(q3)
# print(q2)
# print(q1)
print(f"\nq1 = {q1.evalf()}")
print(f"q2 = {q2.evalf()}")
print(f"q3 = {q3.evalf()}")




#print("\n###########p=f(q)")
# test if valid
pt = [0,0,0]

pt[0] = L*cos(q1) + N*cos(q1+q2)*cos(q3)
pt[1] = L*sin(q1) + N*sin(q1+q2)*cos(q3)
pt[2] = M + N*sin(q3)

pt = Matrix(pt)

print("\nf(q):")

pprint(pt)
print()

input("\n\npress Enter")




#### E3
os.system(CLEAR_STR)

print("\n\n###########  ES3")

q1,q2,q3 = symbols('q1., q2., q3.')   #"." added due to visual bugs in pprint
L, M, N = symbols('L, M, N')

pt[0] = L*cos(q1) + N*cos(q1+q2)*cos(q3)
pt[1] = L*sin(q1) + N*sin(q1+q2)*cos(q3)
pt[2] = M + N*sin(q3)

pt = Matrix(pt)

print("\nf(q):")
pprint(pt)
print()
J = pt.jacobian((q1,q2,q3))

print("\nJf(q):")
pprint(J)

print("\ndet(Jf(q)):")
print(det(J))

print("\neig(Jf([0,pi/2,0])):")
print((J.subs(q1,0).subs(q2,pi/2).subs(q3,0).subs(N,0.5).subs(L,0.5)).eigenvals())

input("\n\npress Enter")




##### E4

os.system(CLEAR_STR)

print("\n\n###########  ES4")

Ji = J.inv()
qs = [-pi/4, pi/4, pi/4]
qg = [0, 0, pi/4]
p_d0 = Matrix([1,-1,0])
#pprint(Ji)
Ji = Ji.subs(q1, qs[0]).subs(q2, qs[1]).subs(q3, qs[2])

print("\nJf(q)^-1")
pprint(Ji)

q_d0 = (Ji@p_d0)
print("\nq_dot(0) = Jf(q(0))^-1*p_dot(0)")
q_d0[1] = q_d0[1].simplify()
pprint(q_d0)
print("\n\n")
# print("\n###########q(t)")


t = Symbol("t")
T = Symbol("T")
T_val = 2
qt = [0,0,0]


vit = [q_d0[0], q_d0[1], q_d0[2]]
for i in range(2):
    Dq = qg[i]-qs[i]
    d = 0
    c = q_d0[i]*T_val/Dq
    b = 3-2*c
    a = 1-c-b
    qt[i] = qs[i] + Dq*(a*(t/T)**3 + b*(t/T)**2+ c*(t/T) + d)
    #s = qs[i] + Dq*(a*(t/T)**3 + b*(t/T)**2+ c*(t/T) + d)
    
    qt[i] = qt[i].subs(T,T_val).subs(N,0.5).subs(L,0.5).subs(M,0.5)
    print(f"\n\nq{i+1}(t):")
    #pprint(qt[i]) # parametric version
    print(qt[i])
    plot(qt[i], (t,0,T_val), title=f"q{i+1}(t)")
    d = diff(qt[i],t,1)
    print(f"\nq{i+1}'(t):")
    print(d)
    plot(d, (t,0,T_val), title=f"q'{i+1}(t)")
    dd = diff(qt[i],t,2)
    print(f"\nq{i+1}''(t):")
    print(dd)
    plot(dd, (t,0,T_val), title=f"q''{i+1}(t)")

Ji = J.inv()
Ji = Ji.subs(q1, qt[0]).subs(q2, qt[1]).subs(q3, qs[2])#qt[2])


print("\n\nJf^-1[3,:]:")
pprint(Ji.subs(t,1)[2,:])

#pprint(J.subs(q3, qs[2]))