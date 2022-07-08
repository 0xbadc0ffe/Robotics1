import numpy as np 
import os
import time
import platform
import copy
import json
from sympy import *
import pyperclip

if platform.system() == 'Windows':
    CLEAR_STR = "cls" 
else:
    CLEAR_STR = "clear"
    
    
print(f"Currently running on: {platform.system()}")


def mat_str(mat, numspace=4 ,trunc=False, large=False, simplification=True, tolerance=1e-10, sc_notation=True):
    # convert matrix in a string, eventually truncating its values (or making them int with
    # trunc = int or trunc = "int". The option large is made to fit any value length ... but it's large ...
    # numspace parameter specify the number of space IN which print the value, or another spacing in large mode

    if simplification:
        mat = nsimplify(mat,tolerance=tolerance,rational=True)

    if sc_notation:
        # s,c act as aliases to sin,cos
        s = sin
        c = cos

        # the intended short-hand functions used to replace sin/cos later on
        s1 = Function('s')
        c1 = Function('c')
        mat = mat.subs({s:s1,c:c1})


    res = ""
    if isinstance(trunc, int):
        if numspace < trunc + 3:
            numspace = trunc + 4

    if large:
        distances=[]
        for row in mat:
            for i in range(len(row)):
                le = len(str(row[i]))
                if len(distances) < len(row):
                    distances.append(le)
                else:
                    if le > distances[i]:
                        distances[i] = le
    for row in mat:
        res += "["
        for i in range(len(row)):
            e = row[i]
            if trunc == "int" or trunc == int:
                e = int(e)
            '''
            elif isinstance(trunc, int) and trunc > 0:
                if e - int(e) >= float('0.'+'9'*trunc):
                    # small correction for x.9999... elements, eg e=x.999314 trunc=3 => e= x+1
                    e = int(e) + 1
                e = truncate(e, trunc)
                # removing the minus if e contains only zeros
                if e.replace('0','') == '-.':
                    e = e[1:]
                # if x.000 => x
                if '.'+'0'*trunc in e:
                    e = e[: -trunc - 1]
            '''      
            if large:
                res += " "*(distances[i]+ numspace - len(str(e))) + f"{e}"
            else:
                res += " "*(numspace - len(str(e))) + f"{e}"
        res += " "*(2) + "]\n"
    return res

def print_mat(mat, numspace=4 ,trunc=False, large=False):
    # stringify the matrix and print it
    print(mat_str(mat, numspace=numspace, trunc=trunc, large=large))

def truncate(f, n):
    #Truncates/pads a float f to n decimal places without rounding
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

def hom_transf_matrix( joint1, joint2 ):
    # compute homogeneous transformation matrix between consecutive joints frame vectors
    # TODO, possible only if implementing joint reference parameters (from which to derive dh parameters)
    return 

def gen_dh_tabel(joint_list):
    # generates the Denavit–Hartenberg parameters table from the joint obj list
    # TABLE ENTRIES: theta, alpha, a, d
    dh_table = []
    for joint in joint_list:
        dh_line = [joint.theta, joint.alpha, joint.a, joint.d]
        dh_table.append(dh_line)
    dh_table = np.array(dh_table)
    return dh_table

def gen_hom_matrix_from_table(index, dh_table):
    # generate homogeneous transformation matrix from dh_table and joint index
    i = index
    hom_mat = np.array([[np.cos(dh_table[i,0]), -np.sin(dh_table[i,0]) * np.cos(dh_table[i,1]), np.sin(dh_table[i,0]) * np.sin(dh_table[i,1]), dh_table[i,2] * np.cos(dh_table[i,0])],
                      [np.sin(dh_table[i,0]), np.cos(dh_table[i,0]) * np.cos(dh_table[i,1]), -np.cos(dh_table[i,0]) * np.sin(dh_table[i,1]), dh_table[i,2] * np.sin(dh_table[i,0])],
                      [0, np.sin(dh_table[i,1]), np.cos(dh_table[i,1]), dh_table[i,3]],
                      [0, 0, 0, 1]])  
    return hom_mat

def input_joint_list(joint_list=[], json_file_name="joints.json", options=None):
    # a simple bash interface used to take the joints parameters in input
    global baseframe, effectorframe, b_e_changed
    print("\nHi, press Enter to start ...\n")
    input()
    os.system(CLEAR_STR)
    #joint_list = []

    if options is None:
        options = { "copy on clip": True,
                    "simplificiation policy": "ending",
                    "simplification type": "trigonometric",
                    "reuse home_b_e when evaluating": False,
                    "print interm": False,
                }

    hom_b_e = None
    ev_hom_b_e = None
    Ja = None
    Ja_v = None

    while(True):
        ans = { "0": "compute", "C": "compute", "c": "compute", "compute": "compute",
                "1": "add", "A": "add", "a": "add", "add": "add",
                "2": "status", "I": "status", "i": "status", "status": "status",
                "3": "edit", "E": "edit", "e": "edit", "edit": "edit",
                "4": "baseframe", "B": "baseframe", "b": "baseframe", "baseframe": "baseframe",
                "5": "effectorframe", "F": "effectorframe", "f": "effectorframe", "effectorframe": "effectorframe",
                "6": "save", "S": "save", "s": "save", "save": "save",
                "7": "import", "I":"import", "i":"import", "import":"import",
                "8": "eval", "V": "eval", "v":"eval", "evaluate":"eval",
                "9": "options", "O":"options", "o":"options", "options":"options",
                "10": "close", "X":"close", "x":"close", "close":"close",
                    }

        bool_ans = { "1": True, "Y": True, "y": True, "yes": True,
                     "0": False, "N": False, "n": False, "no": False,
        }
        
        print(f"\nCurrent number of joints: {len(joint_list)}")
        #print("\nWould you like to add a new joint?    \n\n")
        print("\nPlease select between the following operations:   \n\n")
        print("0/C/compute:    compute Transf. from joints      1/A/add:       add joint to joints list\n\n")
        print("2/I/status:     show status                      3/E/edit:      edit/remove joint\n\n")
        print("4/B/baseframe:  change base frame                5/F/effector:  change effector frame\n\n")
        print("6/S/save:       save joints list into file       7/I/import:    import joints list from file\n\n")
        print("8/V/evaluate:   evaluate                         9/O/options:   computation options\n\n")
        print("10/X/close:     close program\n\n\n")
        inp = input()
        try:
            sw = ans[inp]

            if sw == "compute":
                if not len(joint_list):
                    os.system(CLEAR_STR)
                    print("\nEmpty joints list!")
                    input("\n\nPress Enter to return\n\n")
                    os.system(CLEAR_STR)
                    continue
                else:
                    os.system(CLEAR_STR)
                    simpp = options["simplificiation policy"]
                    cpp = options["copy on clip"]
                    hom_b_e, Ja = compute_all(joint_list, trunc=3, large=True, simplification_policy=simpp, clipcopy=cpp)
                    input("\n\nPress Enter to return\n\n")
                    os.system(CLEAR_STR)
                    continue

            elif sw == "add":
                # Addin a joint
                os.system(CLEAR_STR)
                index = len(joint_list) + 1
                print(f"\nJoint n° {index}:\n")
                print("For each entry input or \"var-variable_name\" if unknown\n")
                while True: 
                    try:
                        theta = input(f"\nTheta{index} (degrees): ")
                        if theta.startswith("var-"):
                            theta = theta[4:]
                        else:
                            theta = float(theta)
                        break
                    except ValueError:
                        print("\nWrong format\n")
                while True: 
                    try:  
                        alpha = input(f"\nAlpha{index} (degrees): ")
                        if alpha.startswith("var-"):
                            alpha = alpha[4:]
                        else:
                            alpha = float(alpha)
                        break
                    except ValueError:
                        print("\nWrong format\n")
                while True: 
                    try:
                        a = input(f"\nA{index} (cm): ")
                        if a.startswith("var-"):
                            a = a[4:]
                        else:
                            a = float(a)
                        break
                    except ValueError:
                        print("\nWrong format")
                while True: 
                    try:
                        d = input(f"\nD{index} (cm): ")
                        if d.startswith("var-"):
                            d = d[4:]
                        else:
                            d = float(d)
                        break
                    except ValueError:
                        print("\nWrong format\n")
                    
                while True:
                    inp = input("\n\nConfirm?   1/Y/yes: yes    0/N/no: no\n\n")
                    try:
                        sw = bool_ans[inp]
                        break
                    except KeyError:
                        print("\n\nPlease use only the given possible answers")
                        continue
                if sw:
                    joint = Joint(index, theta, alpha, a, d)
                    joint_list.append(joint)
                
                os.system(CLEAR_STR)
                continue

            elif sw == "status":
                os.system(CLEAR_STR)
                if len(joint_list) > 0:
                    print_joint_list(joint_list)
                else:
                    print("\nJoint list is empty ...")
                print("\n\n[Base Frame transformation]:\n")
                print_mat(baseframe, trunc=3)
                print("\n\n[Effector Frame transformation]:\n")
                print_mat(effectorframe, trunc=3)
                if hom_b_e is not None:
                    print("\n\n[Last computed hom_b_e]:\n")
                    print_mat(hom_b_e, trunc=trunc, large=True)
                if Ja is not None:
                    print("\n\n[Last computed Ja]:\n")
                    print_mat(np.array(Ja), trunc=trunc, large=True)
                if ev_hom_b_e is not None:
                    print("\n\n[Last computed ev_hom_b_e]:\n")
                    print_mat(ev_hom_b_e, trunc=trunc, large=True)
                if Ja_v is not None:
                    print("\n\n[Last computed Ja_v]:\n")
                    print_mat(np.array(Ja_v), trunc=trunc, large=True)


                input("\n\nPress Enter to return\n\n")
                os.system(CLEAR_STR)
                continue
            
            elif sw == "edit":
                os.system(CLEAR_STR)
                if len(joint_list) > 0:
                    print_joint_list(joint_list)
                    while True:
                        op = input ("\n\nChoose operation [e: edit | r: remove | x: exit]: ")
                        if op == "r" or op == "e" or op == "x":
                            break

                    if op == "x":
                        os.system(CLEAR_STR)
                        continue
               
                    while True:
                        if op == "r":
                            n = input("\n\nJoint to remove: ")
                        else:
                            n = input("\n\nJoint to edit: ")
                        try:
                            n = int(n)
                            if n-1 < len(joint_list):
                                break
                            else:
                                print("Joint not in list")
                        except ValueError:
                            print("Wrong Fromat") 
                    
                    if op == "r":                  
                        joint_list.pop(n-1)
                        # renaming joints
                        for i in range(len(joint_list)):
                            joint_list[i].num = i + 1
                        os.system(CLEAR_STR)
                        print("\nJoint successfully removed!")
                        input("\n\nPress Enter to return\n\n")
                        os.system(CLEAR_STR)
                        continue
                    elif op == "e":
                        os.system(CLEAR_STR)
                        print(f"\nJoint n° {n}:\n")
                        print("For each entry input its value or x if unknown\n")
                        while True: 
                            try:
                                theta = input(f"\nTheta{n} (degrees): ")
                                if theta.startswith("var-"):
                                    theta = theta[4:]
                                else:
                                    theta = float(theta)
                                break
                            except ValueError:
                                print("\nWrong format\n")
                        while True: 
                            try:  
                                alpha = input(f"\nAlpha{n} (degrees): ")
                                if alpha.startswith("var-"):
                                    alpha = alpha[4:]
                                else:
                                    alpha = float(alpha)
                                break
                            except ValueError:
                                print("\nWrong format\n")
                        while True: 
                            try:
                                a = input(f"\nA{n} (cm): ")
                                if a.startswith("var-"):
                                    a = a[4:]
                                else:
                                    a = float(a)
                                break
                            except ValueError:
                                print("\nWrong format")
                        while True: 
                            try:
                                d = input(f"\nD{n} (cm): ")
                                if d.startswith("var-"):
                                    d = d[4:]
                                else:
                                    d = float(d)
                                break
                            except ValueError:
                                print("\nWrong format\n")
                            
                        while True:
                            inp = input("\n\nConfirm?   1/Y/yes: yes    0/N/no: no\n\n")
                            try:
                                sw = bool_ans[inp]
                                break
                            except KeyError:
                                print("\n\nPlease use only the given possible answers")
                                continue
                        if sw:
                            joint = Joint(n, theta, alpha, a, d)
                            joint_list[n-1]=joint
                            os.system(CLEAR_STR)
                            print("\nJoint successfully edited!")
                            input("\n\nPress Enter to return\n\n")
                        os.system(CLEAR_STR)
                        continue

                        
                else:
                    print("\nJoint list is empty ...")
                    input("\n\nPress Enter to return\n\n")
                    os.system(CLEAR_STR)
                    continue           
            
            elif sw == "baseframe":
                # Takes hom. transformation matrix from baseframe to starting joint in input
                os.system(CLEAR_STR)
                print("\nDescribe the homogeneous transformation matrix from base frame to the starting joint\n\n")
                frame_mat = []
                for i in range(4):
                    while(True):
                        r_i = input(f"[Row {i+1}]  |  ")
                        ans = input("Enter to confirm, \"r\" to repeat\n")
                        if ans != "r":
                            try:
                                r_i = r_i.split()
                                if len(r_i) != 4 :
                                    print("Enter 4 numbers      e.g 1 2.4 0 0.93\n")
                                    continue
                                frame_mat_line = []
                                for number in r_i:
                                    frame_mat_line.append(float(number))
                                frame_mat.append(frame_mat_line)
                                break
                            except ValueError:
                                print("Wrong number format\n")
                                continue

                baseframe = np.array(frame_mat)
                b_e_changed = True
                os.system(CLEAR_STR)
                print()
                print_mat(baseframe, trunc=3)
                input("\nPress Enter to return ...")
                os.system(CLEAR_STR)
                continue
            
            elif sw == "effectorframe":
                # Takes hom. transformation matrix from ending joint to effector frame in input
                os.system(CLEAR_STR)
                print("\nDescribe the homogeneous transformation matrix from the ending joint to the effector frame\n\n")
                frame_mat = []
                for i in range(4):
                    while(True):
                        r_i = input(f"[Row {i+1}]  |  ")
                        ans = input("Enter to confirm, \"r\" to repeat\n")
                        if ans != "r":
                            try:
                                r_i = r_i.split()
                                if len(r_i) != 4 :
                                    print("Enter 4 numbers      e.g 1 2.4 0 0.93\n")
                                    continue
                                frame_mat_line = []
                                for number in r_i:
                                    frame_mat_line.append(float(number))
                                frame_mat.append(frame_mat_line)
                                break
                            except ValueError:
                                print("Wrong number format\n")
                                continue

                effectorframe = np.array(frame_mat)
                b_e_changed = True
                os.system(CLEAR_STR)
                print()
                print_mat(effectorframe, trunc=3)
                input("\nPress Enter to return ...")
                os.system(CLEAR_STR)
                continue
              
            elif sw == "save":
                data = {}

                if not b_e_changed:
                    data["baseframe"] = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
                    data["effectorframe"] = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
                else:
                    if baseframe is not None :
                        data["baseframe"] = baseframe.tolist()
                    else:
                        data["baseframe"] = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
                    if effectorframe is not None:
                        data["effectorframe"] = effectorframe.tolist()
                    else:
                        data["effectorframe"] = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

                data["joints"] = []
                for j in joint_list:
                    var = {"th":j.theta_v, "al":j.alpha_v, "a":j.a_v, "d":j.d_v}
                    data["joints"].append(var)

                with open(json_file_name, "w+") as jfile:
                    json.dump(data, jfile, indent=4)
                os.system(CLEAR_STR)
                print("\nData succesfully saved!\n")
                input("\nPress Enter to return ...")
                os.system(CLEAR_STR)
                continue

            elif sw == "import":
                with open(json_file_name, "r") as jfile:
                    data = json.load(jfile)
                    effectorframe = np.array(data["effectorframe"])
                    baseframe = np.array(data["baseframe"])
                    joint_list_vars = data["joints"]
                joint_list = []
                for num,j in enumerate(joint_list_vars):
                    joint_list.append(Joint(num, *list(j.values()))) #give_deg2rad=False))
                os.system(CLEAR_STR)
                print()
                print("Data succesfully imported!\n")
                if len(joint_list) > 0:
                    print_joint_list(joint_list)
                else:
                    print("\nJoint list is empty ...")
                print("\n\n[Base Frame transformation]:")
                print_mat(baseframe, trunc=3)
                print("\n\n[Effector Frame transformation]:")
                print_mat(effectorframe, trunc=3)
                input("\nPress Enter to return ...")
                os.system(CLEAR_STR)
                continue

            elif sw == "eval":
                if not len(joint_list):
                    os.system(CLEAR_STR)
                    print("\nEmpty joints list!")
                    input("\n\nPress Enter to return\n\n")
                    os.system(CLEAR_STR)
                    continue
                else:
                    os.system(CLEAR_STR)
                    ev_hom_b_e, Ja_v = compute_all(joint_list, trunc=3, large=True, evaluate=True)
                    input("\n\nPress Enter to return\n\n")
                    os.system(CLEAR_STR)
                    continue

            elif sw == "options":
                os.system(CLEAR_STR)
                print("\nNot yet avaible \(t >t)/ :c \n")
                input("\n\nPress Enter to return\n\n")
                os.system(CLEAR_STR)

            elif sw == "close":
                close()

        except KeyError:
            os.system(CLEAR_STR)
            print("\nPlease use only the given possible answers\n")
            input("\n\nPress Enter to return\n\n")
            os.system(CLEAR_STR)
            continue
            

def close(timesl=1):
    # close the program
    os.system(CLEAR_STR)
    print("\n\n\n           Bye         ,(è >è)/\n\n\n")
    time.sleep(timesl)
    os.system(CLEAR_STR)
    exit()

def print_joint_list(joint_list):
    # print the joint list
    print("\nJoint List:\n")
    for j in joint_list:
        print()
        print(j)
        print()

def compute_all(joint_list, trunc=3, large=False, differentiate=True, evaluate=False, simplification_policy="ending", simplification_type="trigonometric", tolerance=1e-10, clipcopy=True):
    # Firstly it computes all the relative frames hom. transformation [Ai-1->i]
    # Then it generate the final transformation from the starting joint frame to the last joint frame
    # IF a baseframe or an effectorframe are given this function also computes the Hom transformation between 
    # base and effector frames

    # simplification_policy: ending| all| None,                 decide when to compute simplifications
    # simplification_type:   trigonometric | general            type of simplification
    # tollerance                                                for numeric approx (but nsimplify gives problems)
    print_joint_list(joint_list)

    if evaluate:
        joint_list_ = [] 
        for j in joint_list:
            if j.q_i == j.theta:
                qi = j.theta_v
            else:
                qi = j.d_v
            joint_list_.append(Joint(j.num, j.theta_v, j.alpha_v, j.a_v, j.d_v, compute_all=True, type=qi))
        joint_list = joint_list_

    input("\n\nPress Enter to compute all homogeneous transformation matrices\n\n")
    os.system(CLEAR_STR)
    hom_list = []
    for j in joint_list:
        print(f"\nMatrix A{j.num - 1}->{j.num}(q{j.num})\n")
        print_mat(j.hom_mat, trunc=trunc, large=large)
        hom_list.append(j.hom_mat)
        input("\nPress Enter to show next ...\n")   
        os.system(CLEAR_STR)
    
    print(f"\nTransformation Frame {0} -> Frame {len(joint_list)}:\n")
    hom_0_n = hom_list[0]
    tot = len(hom_list)
    #hom_0_n = nsimplify(hom_0_n,tolerance=tolerance,rational=True)
    if simplification_type == "trigonometric":
        simp = trigsimp
    elif simplification_type == "general":
        simp = simplify
    else:
        simp = lambda x: x

    for n, h in enumerate(hom_list[1:]):
        hom_0_n = hom_0_n @ h
        #hom_0_n = nsimplify(hom_0_n,tolerance=tolerance,rational=True)
        if simplification_policy == "all":
            print(f"Computation & simplification:   {truncate((n+1)/tot*100, 2)} %")
            for i in range(len(hom_0_n[:,0])):
                for j in range(len(hom_0_n[0,:])):
                    #print(hom_0_n[i,j])
                    hom_0_n[i,j] = nsimplify(hom_0_n[i,j],tolerance=tolerance,rational=True)
                    hom_0_n[i,j] = simp(hom_0_n[i,j])
                    #print(hom_0_n[i,j], "\n")

    if simplification_policy == "ending":
        print(f"Simplification in act ... may take a while ...")
        tot = len(hom_0_n[:,0])*len(hom_0_n[0,:])
        for i in range(len(hom_0_n[:,0])):
            for j in range(len(hom_0_n[0,:])):
                print(f"Simplification:   {(i*4+j)/tot*100} %")
                #print(hom_0_n[i,j])
                hom_0_n[i,j] = nsimplify(hom_0_n[i,j],tolerance=tolerance,rational=True)
                hom_0_n[i,j] = simp(hom_0_n[i,j])
                #print(hom_0_n[i,j], "\n")


    os.system(CLEAR_STR)
    print(f"\nTransformation Frame {0} -> Frame {len(joint_list)}:\n")

    if evaluate:
        # TODO: instead of this, avoid numeric simplification if evaluate=True
        for i in range(len(hom_0_n[:,0])):
            for j in range(len(hom_0_n[0,:])):
                hom_0_n[i,j] = hom_0_n[i,j].evalf()


    print_mat(hom_0_n, trunc=trunc, large=large)
    
    global baseframe, effectorframe, b_e_changed
    if b_e_changed:
        input("\nPress Enter to show next ...\n")   
        os.system(CLEAR_STR)
    
        print(f"\nTransformation base frame -> effector frame:\n")
        hom_b_e = baseframe @ hom_0_n @ effectorframe
        print_mat(hom_b_e, trunc=trunc, large=large)
    else:
        hom_b_e = hom_0_n

    pe = hom_b_e[:-1,-1]
    Ja = zeros(pe.shape[0], len(joint_list))
    for k in range(pe.shape[0]): 
        for i, j in enumerate(joint_list):
            try:
                Ja[k,i] = pe[k].diff(j.q_i)
            except ValueError:
                #print(j)
                #print(j.q_i)
                print(f"cannot compute {i}-th row of Ja.\nThis may be caused by a wrong setup of the D-H table.")
                pass

    if clipcopy:
        pyperclip.copy(mat_str(hom_0_n, large=True, simplification=(not evaluate), sc_notation=(not evaluate)))
        print("\n\n+ Data on clipboard! +")
    
    if differentiate:
        print("\n\n\nAnalytic Jacobian:\n\n")
        print_mat(np.array(Ja), trunc=trunc, large=large)
        return hom_b_e, Ja
    else:
        return hom_b_e




class Joint():
    # This class define the joint characteristics given by the D-H parameters

    def __init__(self, num, theta, alpha, a, d, give_deg2rad=True, compute_all=False, type="?"):
        self.num = num                            # joint number
        self.q_i = type
        self.theta = Symbol(f"θ{num}")
        self.alpha = Symbol(f"α{num}")
        self.a = Symbol(f"a{num}")
        self.d = Symbol(f"d{num}")
        self.inst_eval_var = []                   # collects variables with common values to evaluate early
        if give_deg2rad:

            self.theta_v = theta
            #self.theta_v = np.deg2rad(theta)      # angle from xi-1 and xi around zi, from degrees
            if compute_all or (not isinstance(theta, str) and theta % 90 == 0):
                if isinstance(theta, str):
                    self.inst_eval_var.append((self.theta, theta))
                else:
                    self.inst_eval_var.append((self.theta, np.deg2rad(theta)))

            self.alpha_v = alpha
            #self.alpha_v = np.deg2rad(alpha)      # angle from zi-1 and zi around xi, from degrees
            if compute_all or (not isinstance(alpha, str) and alpha % 90 == 0):
                if isinstance(alpha, str):
                    self.inst_eval_var.append((self.alpha, alpha))
                else:
                    self.inst_eval_var.append((self.alpha, np.deg2rad(alpha)))

        else:
            self.theta_v = theta                  # angle from zi-1 and zi around xi
            self.alpha_v = alpha                  # angle from zi-1 and zi around xi
            if compute_all or (not isinstance(theta, str) and (abs(theta % np.pi/2) > np.pi/2-0.01 or abs(theta % np.pi/2) < 0.01)):
                self.inst_eval_var.append((self.theta, theta))
            if compute_all or (not isinstance(alpha, str) and (abs(alpha % np.pi/2) > np.pi/2-0.01 or abs(alpha % np.pi/2) < 0.01)):
                self.inst_eval_var.append((self.alpha, alpha))
        self.a_v = a                              # distance of origin of frame i-1 to origin of frame i along xi-1
        self.d_v = d                              # distance of origin of frame i-1 to origin of frame i along zi-1
        if compute_all or (not isinstance(a, str) and a == 0):
            self.inst_eval_var.append((self.a, a))
        if compute_all or (not isinstance(d, str) and d == 0):
            self.inst_eval_var.append((self.d, d))
        self.hom_mat = None                       # homogeneous transformation matrix from frame i-1 to frame i
        self.gen_hom_matrix()

        if self.q_i == "?":
            if isinstance(theta, str) ^ isinstance(d, str):
                if isinstance(theta, str):
                    self.q_i = self.theta
                if isinstance(d, str):
                    self.q_i = self.d 
            else:
                os.system(CLEAR_STR)
                print("\nAmbiguity for:\n")
                print(self)
                while True:
                    ans = input("\n\nType of joint [p: primatic | r: revolute]: ")
                    if ans == "p":
                        self.q_i = self.d 
                        break
                    elif ans == "r":
                        self.q_i = self.theta
                        break
                    else:
                        print("\nWrong Format")





    def gen_hom_matrix(self):
        # generate the homogeneous transformation matrix
        '''
        hom_mat = np.array([[np.cos(self.theta), -np.sin(self.theta) * np.cos(self.alpha), np.sin(self.theta) * np.sin(self.alpha), self.a * np.cos(self.theta) ],
                      [np.sin(self.theta), np.cos(self.theta) * np.cos(self.alpha), -np.cos(self.theta) * np.sin(self.alpha), self.a * np.sin(self.theta) ],
                      [0, np.sin(self.alpha), np.cos(self.alpha), self.d ],
                      [0, 0, 0, 1]]) 
        '''
        hom_mat = np.array([[cos(self.theta), -sin(self.theta) * cos(self.alpha), sin(self.theta) * sin(self.alpha), self.a * cos(self.theta) ],
                      [sin(self.theta), cos(self.theta) * cos(self.alpha), -cos(self.theta) * sin(self.alpha), self.a * sin(self.theta) ],
                      [0, sin(self.alpha), cos(self.alpha), self.d ],
                      [0, 0, 0, 1]])  

        for i in range(len(hom_mat)):
            for j in range(len(hom_mat[0])): 
                try:

                    hom_mat[i,j] = hom_mat[i,j].subs(self.inst_eval_var).evalf()

                except AttributeError:
                    pass
                    #hom_mat[i,j]
                except SympifyError: 
                    pass      

        self.hom_mat = hom_mat
        return hom_mat

    def __str__(self):
        # string representation
        if isinstance(self.theta_v,float):
            th = truncate(self.theta_v, 4)
        else:
            th = self.theta_v
        if isinstance(self.alpha_v,float):
            al = truncate(self.alpha_v, 4)
        else:
            al = self.alpha_v
        if isinstance(self.a_v,float):
            a = truncate(self.a_v, 4)
        else:
            a = self.a_v
        if isinstance(self.d_v,float):
            d = truncate(self.d_v, 4)
        else:
            d = self.d_v
        #th = truncate(self.theta, 4)
        #al = truncate(self.alpha, 4)
        #a = truncate(self.a, 1)
        #d = truncate(self.d, 1)
        res = f"[Joint {self.num}]  Theta: {th} [rad]   Alpha: {al} [rad]   A: {a} [cm]   D: {d} [cm]"

        res = f"[Joint {self.num}]  Theta: {th} [rad]   Alpha: {al} [deg]   A: {a} [cm]   D: {d} [cm]     q{self.num}:{self.q_i}"
        return res




'''
Another way to compute stuff

# create a Joint obj
joint1 = Joint(1, 90, 90, 0, 2)

# print a matrix truncating at 3rd decimal
print_mat(joint1.hom_mat, trunc=3)

# adding joint to joint list
joint_list = [joint1]

# generate dh table
dh_table = gen_dh_tabel(joint_list)
print(dh_table)

# generate and print the homogeneous transf. matrix from the dh table and the hom_mat field of the Joint objects
print_mat(gen_hom_matrix_from_table(0, dh_table), trunc=3)            
'''



###### MAIN

if __name__ == "__main__":
    
    # Initializing Base and Effector frames
    baseframe = np.identity(4)
    effectorframe = np.identity(4)
    # b_e_changed keeps track of baseframe or effectorframe changes
    b_e_changed = False 

    # generate joint list and data
    joint_list = input_joint_list()


    '''
    # compute all transformations
    hom_0_n = compute_all(joint_list, trunc=3, large=True, simplification_policy="ending")

    input("\nPress Enter to evaluate\n\n")
    os.system(CLEAR_STR)

    # evaluate
    hom_0_n = compute_all(joint_list, trunc=3, large=True, evaluate=True)
    '''
    
    input("\nPress Enter to exit\n\n")
    close()


    
