from sympy import *
import os

os.system("cls")

t = Symbol("t")
T = Symbol("T")

q0 = Symbol("q0")
qf = Symbol("qf")
v0 = Symbol("v0")
vf = Symbol("vf")
a0 = Symbol("a0")
af = Symbol("af")

# General Quintic polynomial
#s = 6*(t/T)**5 -15*(t/T)**4 + 10*(t/T)**3  # goes from 0 to 1 in T secs

s = (1-(t/T))**3*(q0 + (3*q0 + v0*T)*(t/T) + 0.5*(a0*T**2 + 6*v0*T + 12*q0)*(t/T)**2) + (t/T)**3*(qf + (3*qf - vf*T)*(1-(t/T)) + 0.5*(af*T**2 - 6*vf*T + 12*qf)*(1-(t/T))**2)


# parameters
q0_val = -60
qf_val = 30
v0_val = 45
vf_val = 0
a0_val = 0
af_val = 0
T_val = 2

s = s.subs(q0,q0_val)
s = s.subs(qf,qf_val)
s = s.subs(v0,v0_val)
s = s.subs(vf,vf_val)
s = s.subs(a0,a0_val)
s = s.subs(af,af_val)


max_d = 3 # maximum derivative plotted (maxima and minima are computed until max_d -1)

# ' symbols
sy = []
for i in range(max_d+1):
    sy.append("\'"*i)

d = s
for k in range(max_d+1):
    s_prev = d
    d = diff(s,t,k)
    z = solve(d,t)
    if k > 0:
        for zi in z:
            print(f"\nmax-min:  t = {zi} =>  ", s_prev.subs(t,zi),"  =>  ", s_prev.subs(t,zi).evalf(), "  =>  ", s_prev.subs(t,zi).subs(T,T_val).evalf(),"\n")
    print(f"\n\n\n{k})  ++++++++++")
    print(f"\ns{sy[k]}(t) = ", d)
    print("\nzeros: ", z)
    plot(d.subs(T,T_val), (t,0,T_val), title=f"s{sy[k]}(t)")


print("\n\n")