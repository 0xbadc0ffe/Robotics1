from sympy import *
import numpy as np

x, theta = symbols("x, θ")



#### Compute Rif from Ri and Rf
Ri = [
    [1,  0,  0],
    [0,  1,  0],
    [0,  0,  1]
]

Ri = Matrix(Ri)

Rf = [
    [0,  1,  0],
    [1,  0,  0],
    [0,  0,  1]
]

Rf = Matrix(Rf)

R = Ri.T@Rf    #Rif


#### if Rif Given

Ri = [
    [1,  0,  0],
    [0,  1,  0],
    [0,  0,  1]
]

R = [

    [0,          -1,  0        ],
    [-1/sqrt(2),  0,  1/sqrt(2)],
    [-1/sqrt(2),  0, -1/sqrt(2)]

]

R = Matrix(R)

print("\nRif:")
pprint(R)



# Compute theta_f and r

th = acos( (R[0,0] + R[1,1] + R[2,2] - 1)/2 ).evalf(9)

rm = [
    [R[2,1] - R[1,2]],
    [R[0,2] - R[2,0]],
    [R[1,0] - R[0,1]]
]
rm = Matrix(rm)

r = 1/(2*sin(theta))*rm

print("\nr:")
pprint(r)
print("\ntheta_f:")
print(th)
print("\nr(theta_f):")
r = r.subs(theta, th).evalf(9)
pprint(r)
print()


# Compute R^i(t)

Sr =[
    [0, -r[2], r[1]],
    [r[2], 0, -r[0]],
    [-r[1], r[0], 0],
]

Sr = Matrix(Sr)


Rit = r @ r.T + (Identity(3) - r @ r.T)*cos(theta) + Sr*sin(theta)

# Re = Ri*R^i(t)
Re = Ri@Rit
decs = 4
Rei = np.trunc(np.array(Re.subs(theta, 0))*10.0**decs)/(10.0**decs)
Ref = np.trunc(np.array(Re.subs(theta, th))*10.0**decs)/(10.0**decs)
print("\nRe(0):")
print(Rei)
print("\nRe(theta_f):")
print(Ref)


#Rg = r @ r.T + (Identity(3) - r @ r.T)*cos(theta) + Sr*sin(theta)
#Rg = Rg.subs(theta, th)
#print(nsimplify(np.array(Rg), tolerance=1e-5, rational=True))
##Rn = np.trunc(np.array(Rg)*10.0**decs)/(10.0**decs)
#print(Rn)
#print(np.around(np.array(Rg),3))
#input()



