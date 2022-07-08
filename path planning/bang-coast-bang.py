from sympy import *
import os

os.system("cls")
t = Symbol("t")

# parameters
q0 = 0
qf = 3
v0 = -2             # eg exam September 11, 2020, ex 5
vmax = "inf" #2     # can be a number or "inf" (symmetric bang-bang)
amax = 4            # NB: if fmax is given => amax = fmax/mass

L = qf - q0   # abs(qf - q0) ?
if vmax != "inf" and (L <= vmax**2/amax or v0 > vmax):  # abs(v0) > vmax?
    if L <= vmax**2/amax:
        print(f"\nL <= vmax**2/amax:   {L} <= {vmax**2/amax}")
    print("\nConditions not satisfied")
    exit()
if vmax != "inf":
    print(f"\nL > vmax**2/amax:   {L} > {vmax**2/amax}\n")

# Slide version (simpler)
#T = (L*amax + vmax**2)/(amax*vmax)    # = L/vmax + vmax/amax
#s_lift = amax * t**2/2
#s_coast = vmax*t - (vmax)**2/(2*amax) 
#s_land = -amax*(t-T)**2/2 +vmax*T - (vmax)**2/amax

if vmax != "inf":
    Tsi = (vmax-v0)/amax
    Ts = (vmax)/amax
    T = L/vmax + Tsi*(vmax-v0)/(2*vmax) +Ts/2

    print(f"\nTsi: {Tsi}")
    print(f"\nTs: {Ts}")
    print(f"\nT: {T}")

    s_lift = amax * t**2/2 + v0*t
    #s_end = amax*Tsi**2/2 + v0*Tsi 

    s_coast = vmax*t - (vmax-v0)**2/(2*amax)
    #s_end = vmax*(T-Tsi) - (vmax-v0)**2/(2*amax)

    s_land = -amax*(t-T)**2/2 +vmax*T - (vmax)**2/(2*amax) - (vmax-v0)**2/(2*amax)

    s = Piecewise((s_lift,t<=Tsi),(s_coast, And(Tsi < t, t<= T-Ts) ),(s_land, And(T-Ts < t , t <= T)))
    #plot(s,(t,0,T), title="q(t)")

else: # symmetric bang-bang
    L = L - amax*(v0/amax)**2/2 + v0*v0/amax    # L = x(T)-x(Tz), so we have to add de displacement 
                                                # x(Tz) = s_lift(t=-v0/amax) to our previous L
    Ts = sqrt(L/amax) - v0/amax
    T = 2*sqrt(L/amax) - v0/amax
    print(f"\nTz: {-v0/amax}")
    print(f"\nTs: {Ts}")
    print(f"\nT: {T}")

    s_lift = amax*t**2/2 + v0*t
    vmax = amax*Ts + v0
    s_land = -amax*(t-T)**2/2 +vmax*T - (vmax)**2/(2*amax) - (vmax-v0)**2/(2*amax)
    s = Piecewise((s_lift,t<=Ts),(s_land, And(Ts < t , t <= T)))

    print(f"\nVmax =  {vmax}")


print(f"\ns(T) = {s.subs(t, T)}\n")

max_d = 2 # maximum derivative plotted (maxima and minima are computed until max_d -1)

# ' symbols
sy = []
for i in range(max_d+1):
    sy.append("\'"*i)

d = s
for k in range(max_d+1):
    s_prev = d
    d = diff(s,t,k)
    try:
        z_chk = False
        z = solve(d,t)
        z_chk = True
    except NotImplementedError:
        pass
    if k > 0 and z_chk:
        for zi in z:
            print(f"\nmax-min:  t = {zi} =>  ", s_prev.subs(t,zi),"  =>  ", s_prev.subs(t,zi).evalf(), "  =>  ", s_prev.subs(t,zi).evalf())
    print(f"\n\n\n{k})  ++++++++++")
    print(f"\ns{sy[k]}(t) = ", d)
    if z_chk:
        print("\nzeros: ", z)
    plot(d, (t,0,T), title=f"s{sy[k]}(t)")


print("\n\n")