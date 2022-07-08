import numpy as np
import os
import platform
import matplotlib.pyplot as plt

if platform.system() == 'Windows':
    CLEAR_STR = "cls" 
else:
    CLEAR_STR = "clear"


def cos(x, deg=False):
    if deg:
        x = np.deg2rad(x)
    return np.cos(x)

def sin(x, deg=False):
    if deg:
        x = np.deg2rad(x)
    return np.sin(x)

def tan(x, deg=False):
    if deg:
        x = np.deg2rad(x)
    return np.tan(x)

def inv( mat):
    mat = np.matrix(mat)
    return mat.I


# Direct kinematics function
# fr examples:  
# fr  = np.array([ L1*cos(q[0]) + L2*cos(q[1])+ L3*cos(q[2]), L1*sin(q[0]) + L2*sin(q[1]) + L3*sin(q[2]), q[2]-q[1]])
def fr(q):
    #fr  = np.array([ L1*cos(q[0]) - q[1]*sin(q[0])+ L2*cos(q[0]-q[2]+np.pi/2) + L3*sin(q[0]-q[2]+np.pi/2), L1*sin(q[0]) - q[1]*cos(q[0])+ L2*sin(q[0]-q[2]+np.pi/2) - L3*cos(q[0]-q[2]+np.pi/2), q[0]-q[2]+np.pi/2])
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    N=0.5
    M=0.5
    L=0.5
    fr  = np.array([L*cos(q1) + N*cos(q1+q2)*cos(q3), L*sin(q1) + N*sin(q1+q2)*cos(q3), M + N*sin(q3)])
    fr = fr[np.newaxis].T
    return fr

# Jacobian Matrix
# Jr examples:
# Jr = np.matrix([[ -L1*sin(q[0]), -L2*sin(q[1]), -L3*sin(q[2]) ], [ L1*cos(q[0]), L2*cos(q[1]), L3*cos(q[2]) ], [0, -1, 1]])
def J(q):
    q = np.squeeze(np.asarray(q))  
    q1 = q[0]
    q2 = q[1]
    q3 = q[2]
    N=0.5
    M=0.5
    L=0.5                                                                                                                                                 
    # Jacobian of fr
    #Jr = np.matrix([[ -L1*sin(q[0]) -q[1]*cos(q[0]) - L2*sin(q[0]-q[2]+np.pi/2) + L3*cos(q[0]-q[2]+np.pi/2), -sin(q[0]),  L2*sin(q[0]-q[2]+np.pi/2) - L3*cos(q[0]-q[2]+np.pi/2) ], [ L1*cos(q[0]) +q[1]*sin(q[0]) +L2*cos(q[0]-q[2]+np.pi/2) + L3*sin(q[0]-q[2]+np.pi/2), -cos(q[0]),  -L2*cos(q[0]-q[2]+np.pi/2) - L3*sin(q[0]-q[2]+np.pi/2)  ], [1, 0, -1]])
    Jr = np.matrix([
        [-L*sin(q1)-N*sin(q1+q2)*cos(q3), -N*sin(q1+q2)*cos(q3), -N*sin(q3)*cos(q1+q2) ],
        [-L*cos(q1)+N*cos(q1+q2)*cos(q3),  N*cos(q1+q2)*cos(q3), -N*sin(q3)*sin(q1+q2) ],
        [0, 0, N*cos(q3)]
    ])
   
    
    return Jr

# q(k+1) = q(k) + Jf(q(k))^-1*[rd -fr(q(k))]
def Qnext(q):
    delta = J(q).I @ [rd -fr(q) ]
    qnext = q + delta.T
    #qnext = np.squeeze(np.asarray(qnext))
    return qnext

# q(k+1) = q(k) + alpha * Jf(q(k))^T*[rd -fr(q(k))]
def Qnext_g(q, alpha=0.7):
    delta = alpha* J(q).T @ [rd -fr(q) ]
    qnext = q + delta.T
    #qnext = np.squeeze(np.asarray(qnext))
    return qnext


# Np settings
np.set_printoptions(precision=4, suppress=True)


# Define variables
# intial guess
#q = np.array([np.pi/2, 1, np.pi/2])
q = np.array([-np.pi/4, np.pi/4, np.pi/4])
q = np.array([0.1, np.pi/2, 0])
q = q[np.newaxis].T

#v = { 'q1':0, 'q2': np.pi/2, "q3": np.pi/2}
# eg v['q1']

# Define constants and constraints
L1 = 1
L2 = 1
L3 = 1
precision = 0.001
alpha = 0.97
method = "NR" # GD or NR
step_print = False

# target
#rd = np.array([2, 1, -np.pi/6])
#rd = np.array([0.3, -0.3, 0.7])
rd = np.array([0.5, 0.5, 0.5])
rd = rd[np.newaxis].T


e_track = []
t_track = []
q_track = {"q1":[q[0]],"q2":[q[1]],"q3":[q[2]]}

os.system(CLEAR_STR)
print("\n######### Starting ...")
print("\nq0: ")
print(q)
print("\nrd: ")
print(rd)
print(f"\n\nprecision: {precision}")
e_norm = precision +1
count = 1
print("\n\n\nEnter to start computation")
while e_norm > precision:    
    if step_print:
        input("\n\nShow next")
    os.system(CLEAR_STR)
    print(f"\n######### Loop {count}:\n")
    print(f"\nq: ")
    print(q%np.pi)
    print("\nfr(q): ")
    print(fr(q))
    print("\ne(q): ")
    e = rd -fr(q)
    print(e)
    e_norm = np.linalg.norm(e)
    print(f"\n||e||: {e_norm}")
    e_track.append(e_norm)
    t_track.append(count)
    print("\nJr(q): ")
    Jr = J(q)
    print(Jr)

    if method == "NR":
        # Newthon - Raphson method
        print("\nDeltaq: ")
        print(np.transpose(Jr.I @ [ rd -fr(q) ]))
        q = Qnext(q)

    elif method == "GD":
        # Gradient method
        print("\nDeltaq")
        print(np.transpose(alpha* J(q).T @ [rd -fr(q) ]))
        print(f"\nq{count}: ")
        q = Qnext_g(q, alpha)
        
    print(q)
    q_track["q1"].append(q[0,0])
    q_track["q2"].append(q[1,0])
    q_track["q3"].append(q[2,0])
    count += 1
    #input()

print("\n\nAccuracy achieved\n")

plt.plot(e_track)   #, t_track)
plt.figure()
plt.plot(q_track["q1"])
plt.plot(q_track["q2"])
plt.plot(q_track["q3"])
plt.show()
input()

