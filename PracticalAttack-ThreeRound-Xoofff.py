# write down this file on 10 June 2020
# this file is to verify the result of attacks on 3-round xoofff



#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random
import re


#define parameters

Xsize = 4
Ysize = 3
Zsize = 4
blocksize = 12 * Zsize
rd_xoodoo = 3 # total rounds of xoodoo used in xoofff
n_out = 372 # the number of output blocks of xoofff 371 for z = 4， 757 for z =8
if Zsize == 4:
    var_theta1 = 1
    var_theta2 = 2
    var_rhowest = 3
    var_rhoeast = 2
    var_rolling1 = 1
    var_rolling2 = 3
elif Zsize  == 8:
    var_theta1 = 3
    var_theta2 = 6
    var_rhowest = 5
    var_rhoeast = 4
    var_rolling1 = 3
    var_rolling2 = 5
elif Zsize == 32:
    var_theta1 = 5
    var_theta2 = 14
    var_rhowest = 11
    var_rhoeast = 8
    var_rolling1 = 5
    var_rolling2 = 13


def tryint(s):                       
    try:
        return int(s)
    except ValueError:
        return s

def str2int(v_str):                
    return [tryint(sub_str) for sub_str in re.split('([0-9]+)', v_str)]

def sort_humanly(v_list):    
    return sorted(v_list, key=str2int)



# C = A * B
def multiply_binary_matrix(A,B):
    rA = len(A)
    lA = len(A[0])
    rB = len(B)
    lB = len(B[0])
    if rA != lA:
        return 'error'
    if rB != lB:
        return 'error'
    if rA != rB:
        return 'error'
    C = [[0 for x in range(rA)] for y in range(lA)]
    for rw in range(rA):
        for cl in range(rA):
            for j in range(rA):
                C[rw][cl] ^= A[rw][j]*B[j][cl]
    return C


        
def solve_equations(B):
    # A is augmented matrix, the last column is constants
    row = len(B)
    col = len(B[0])
    print('row = ' + str(row))
    print('col = ' + str(col))
    if row == 1:
        for i in range(col):
            if B[row-1][i] == 1:
                return 1
        return 0
    #col-1 * col-1 matrix 下三角矩阵， 变量个数为col-1
    if row < (col-1):
        #matrix is row * col
        total = row
    else:
        total = col - 1
    for r in range(total):                                        
        if B[r][r] == 0:
            #exchange two row, if all B[.][r] = 0, go on 
            for rw in range(r+1, row):
                if B[rw][r] == 1:
                    tem = [0 for i in range(col)]
                    for cl in range(r, col):
                        tem[cl] = B[r][cl]
                        B[r][cl] = B[rw][cl]
                        B[rw][cl] = tem[cl]
                    break                    
                
        #make all B[rw][r] = 0 rw = r+1,.....
        for rw in range(r+1,row):
            if B[rw][r] == 1:
                for cl in range(r, col):
                    B[rw][cl] = B[r][cl] ^ B[rw][cl]

        #if the last row, solve the equations
        if r == (total-1):
            Guess_K = [0 for i in range(blocksize)]
            
            for i in range(blocksize):
                if B[col-2-i][col-2-i] == 0:
                    print('k'+str(col-2-i) +' is free variable')
                    return 'failed'
                if i == 0:
                    Guess_K[blocksize - 1 - i] = B[col-2][col - 1]
                else:
                    Guess_K[blocksize - 1 - i] = B[col-2-i][col - 1]
                    for j in range(i):
                        Guess_K[blocksize - 1 - i] ^= B[col-2-i][col-2-j] * Guess_K[blocksize - 1 - j]
            print(Guess_K)            
            return Guess_K

def inverse_binary_matrix(A):
    row = len(A)
    col = len(A[0])
    if row != col:
        return "error"
    
    #first inv_A is set as identified matrix
    inv_A = [[0 for y in range(col)] for x in range(row)]
    for rw in range(row):
        inv_A[rw][rw] = 1
    
    B = [[0 for y in range(col)] for x in range(row)]
    for x in range(row):
        for y in range(col):
            B[x][y] = A[x][y]
            
    #step 1: transfer A to a upper triangular matrix, inv_A have the same operations
    for r in range(row-1):
        if B[r][r] == 0:
            for rw in range(r+1, row):
                if B[rw][r] == 1:
                    tem = 0                    
                    for cl in range(col):
                        tem = B[r][cl]
                        B[r][cl] = B[rw][cl]
                        B[rw][cl] = tem
                        
                        tem = inv_A[r][cl]
                        inv_A[r][cl] = inv_A[rw][cl]
                        inv_A[rw][cl] = tem
                    break
        
        for rw in range(r+1, row):
            if B[rw][r] == 1:
                for cl in range(col):
                    B[rw][cl] = B[rw][cl] ^ B[r][cl]
                    inv_A[rw][cl] = inv_A[rw][cl] ^ inv_A[r][cl]

    #step 2: transfer A into identified matrix, inv_A have the same operations
    for cl in range(col-1, 0, -1):
        for rw in range(cl):
            if B[rw][cl] == 1:
                for j in range(col):
                    B[rw][j] = B[rw][j] ^ B[cl][j]
                    inv_A[rw][j] = inv_A[rw][j] ^ inv_A[cl][j]
        
    return inv_A    
         
                            
    
#theta function in xoodoo
def theta_function(A):
    
    # P = A0 + A1 + A2
    P = [[0 for z in range(Zsize)] for x in range(Xsize)]
    for x in range(Xsize):
        for z in range(Zsize):
            P[x][z] = A[x][0][z] ^ A[x][1][z] ^ A[x][2][z]
            
    # E = P <<< (1,3) + P <<< (1,6) for z = 8
    # E = P <<< (1,1) + P <<< (1,2) for z = 4
    E = [[0 for z in range(Zsize)] for x in range(Xsize)]
    for x in range(Xsize):
        for z in range(Zsize):
            E[x][z] = P[(x-1) % Xsize][(z-var_theta1) % Zsize] ^ P[(x-1) % Xsize][(z-var_theta2) % Zsize]

    #A_y = A_y + E
    for x in range(Xsize):
        for z in range(Zsize):
            A[x][0][z] = A[x][0][z] ^ E[x][z]
            A[x][1][z] = A[x][1][z] ^ E[x][z]
            A[x][2][z] = A[x][2][z] ^ E[x][z]
    return A

# rho_{west} function in xoosdoo
def rho_west_function(A):

    B = [[0 for z in range(Zsize)] for x in range(Xsize)]
    C = [[0 for z in range(Zsize)] for x in range(Xsize)]

    for x in range(Xsize):
        for z in range(Zsize):
            #A_1 = A_1 <<< (1,0) 
            B[x][z] = A[(x-1) % Xsize][1][z]
            #A_2 = A_2 <<< (0,5) for z = 8
            #A_2 = A_2 <<< (0,3) for z = 4
            C[x][z] = A[x][2][(z-var_rhowest) % Zsize]

    for x in range(Xsize):
        for z in range(Zsize):
            A[x][1][z] = B[x][z]
            A[x][2][z] = C[x][z]
    return A
            
# iota function in xoodoo
def iota_function(A, const):
    
    for z in range(Zsize):
        A[0][0][z] = A[0][0][z] ^ ((const >> z) & 1)
    return A

# chi function in xoodoo
def chi_function(A):
    
    # B_i = ~A_{i+1}*A_{i+2} ^ A_{i}
    B = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                B[x][y][z] = ((1 ^ A[x][(y+1)%Ysize][z]) & A[x][(y+2)%Ysize][z]) ^ A[x][y][z]

    return B

# rho_{east} function in xoodoo
def rho_east_function(A):
    
    B = [[0 for z in range(Zsize)] for x in range(Xsize)]
    C = [[0 for z in range(Zsize)] for x in range(Xsize)]

    for x in range(Xsize):
        for z in range(Zsize):
            #A_1 = A_1 <<< (0,1)
            B[x][z] = A[x][1][(z-1) % Zsize]
            #A_2 = A_2 <<< (2,4) for z = 8
            #A_2 = A_2 <<< (2,2) for z = 4
            C[x][z] = A[(x-2) % Xsize][2][(z-var_rhoeast) % Zsize]

    for x in range(Xsize):
        for z in range(Zsize):
            A[x][1][z] = B[x][z]
            A[x][2][z] = C[x][z]
    return A

# NLFSR, roll_e
def roll_e_NLFSR_function(A):
    
    #A[0,0] = A[0,1]* A[0,2] + A[0,0]<<<3 + A[0,1]<<<5 + 0x07 for z = 8
    #A[0,0] = A[0,1]* A[0,2] + A[0,0]<<<1 + A[0,1]<<<3 + 0x3 for z = 4
    
    C = [0 for z in range(Zsize)]
    for z in range(Zsize):
        C[z] = (A[0][1][z] & A[0][2][z]) ^ A[0][0][(z-var_rolling1) % Zsize] ^ A[0][1][(z-var_rolling2) % Zsize] ^ ((0x00000007 >> z) & 1)

    for z in range(Zsize):
        A[0][0][z] = C[z]
    # B = A_0 <<< (3,0)
    B = [[0 for z in range(Zsize)] for x in range(Xsize)]
    for x in range(Xsize):
        for z in range(Zsize):
            B[x][z] = A[(x-3) % Xsize][0][z]
    # A_0 = A_1, A_1 = A_2, A_2 = B
    for x in range(Xsize):
        for z in range(Zsize):
            A[x][0][z] = A[x][1][z]
            A[x][1][z] = A[x][2][z]
            A[x][2][z] = B[x][z]
    return A
 
def Xoodoo_function (S, K):
    const = [0x00000060, 0x0000002c, 0x00000380, 0x000000f0, 0x000001a0, 0x00000012]
    A = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                A[x][y][z] = S[x][y][z]
    
    for r in range(rd_xoodoo):
        A = theta_function(A)
        A = rho_west_function(A)
        A = iota_function(A,const[r])
        A = chi_function(A)
        if r < (rd_xoodoo - 1):
            A = rho_east_function(A)   # at the last round, we don't care this operation

    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                A[x][y][z] =  A[x][y][z] ^ K[x][y][z]
    return A

        

def three_round_attack(inv_Ma):
    
    S = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    C = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    K = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    # random choose S
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                S[x][y][z] = random.randint(0,1)

    true_key = [0 for i in range(blocksize)]
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                K[x][y][z] = random.randint(0,1)
                true_key[3 * x + y + 12 * z] = K[x][y][z]

    print("True key:")
    print(true_key)           
                    
    #date stream is SS, output stream of xoofff is CC
    SS = [S for i in range(n_out)]
    CC = [C for i in range(n_out)]
    
    #o = open("test.txt", 'w')
    for n in range(n_out):
        CC[n] = Xoodoo_function(SS[n], K)
        if n < (n_out -1):
            SS[n+1] = roll_e_NLFSR_function(SS[n])

    #list all variables
    Total_vars = []
    #monomials with degree 1, totally 96 items for z = 8
    for i in range (blocksize):
        Total_vars.append([i])
    #monomials with degree 2 totally 96 items for z = 8
    for i in range(blocksize):
        if (i % 3) == 2:
            Total_vars.append([i-2, i])
        else:
            Total_vars.append([i, i+1])
    for i in range(0,blocksize, 3):
        for j in range(i+3, blocksize):
            Total_vars.append([i, j])
            Total_vars.append([i+1, j])
            Total_vars.append([i+2, j])
        
    #monomials with degree 3  totally 8960 items for z = 8
    for i in range(blocksize):
        for j in range(blocksize):
            if (i % 3) == 0:
                if j != i and j != (i + 2):
                    item = Total_vars[i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
            elif (i % 3) == 1:
                if j != i and j != (i - 1) and j != (i + 1):
                    item = Total_vars[i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
            else:                
                if j != i and j != (i - 1) and j != (i - 2):
                    item = Total_vars[i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
    #monomials with degree 4 totally 4464 items for z = 8
    for i in range(blocksize-3):
        for j in range(i - int(i%3) + 3, blocksize):
            if (i % 3) == 0:
                if j != i and j != (i + 1) and j != (i + 2):
                    item = Total_vars[blocksize+i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
            elif (i % 3) == 1:
                if j != i and j != (i - 1) and j != (i + 1):
                    item = Total_vars[blocksize+i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
            else:                
                if j != i and j != (i - 1) and j != (i - 2):
                    item = Total_vars[blocksize+i] + Total_vars[blocksize+j]
                    #item = sort_humanly(item)
                    item.sort()
                    Total_vars.append(item)
    Total_vars = Total_vars[blocksize:] + Total_vars[0:blocksize]
    vars_ln = len(Total_vars)       

    #the state berfore the last chi
    M_chi_r3 = [[]for i in range(blocksize)]
    for i in range(blocksize):
        y = (i%12)%3
        term = [i-y+(i+1)%3,i-y+(i+2)%3]
        term.sort()
        M_chi_r3[i] = [term, [i],[i-y+(i+2)%3]]   # variable
        #M_chi_r3[i].sort()    
    #the state before last iota, note the third round xor 0x00000380

    for z in range(Zsize):
        con = 0x00000380
        if ((con >> z) & 1) == 1:
            M_chi_r3[12*z].append(['1'])   
            
    # the state before the second rho_east
    M_roest_r2 = [[] for i in range(blocksize)]
    for i in range(blocksize):
        for col in range(blocksize):
            if inv_Ma[i][col] == 1:
                for term in M_chi_r3[col]:
                    if term not in M_roest_r2[i]:
                        M_roest_r2[i].append(term)
                    else:
                        M_roest_r2[i].remove(term)    

    #the state before second chi
    M_chi_r2 = [[] for i in range(blocksize)]
    for i in range(blocksize):
        # x0 + x2 + x1*x2
        # put all terms in x0 into new state first

        for term in M_roest_r2[i]:
            if term not in M_chi_r2[i]:
                M_chi_r2[i].append(term)
            else:
                M_chi_r2[i].remove(term)

        y = (i%12)%3
        if y == 0:
            x1 = i + 1
            x2 = i + 2
        if y == 1:
            x1 = i + 1
            x2 = i - 1
        if y == 2:
            x1 = i - 2
            x2 = i - 1
        
        #put all terms in x2 into new state
            
        for term in M_roest_r2[x2]:
            if term not in M_chi_r2[i]:
                M_chi_r2[i].append(term)
            else:
                M_chi_r2[i].remove(term)
        # handle with x1*x2
        
        for tx1 in M_roest_r2[x1]:                
            for tx2 in M_roest_r2[x2]:
                if tx1 != ['1'] and tx2 != ['1']:
                    t = []
                    for term in tx1:
                        if term not in t:
                            t.append(term)
                    for term in tx2:
                        if term not in t:
                            t.append(term)
                    #t = sort_humanly(t)
                    t.sort()
                    if t not in M_chi_r2[i]:
                        M_chi_r2[i].append(t)
                    else:
                        M_chi_r2[i].remove(t)
                        
                if tx1 == ['1'] and tx2 != ['1']:                        
                    if tx2 not in M_chi_r2[i]:
                        M_chi_r2[i].append(tx2)
                    else:
                        M_chi_r2[i].remove(tx2)                        
        

                if tx1 != ['1'] and tx2 == ['1']:                        
                    if tx1 not in M_chi_r2[i]:
                        M_chi_r2[i].append(tx1)
                    else:
                        M_chi_r2[i].remove(tx1)
                
                if tx1 == ['1'] and tx2 == ['1']:                        
                    if tx1 not in M_chi_r2[i]:
                        M_chi_r2[i].append(tx1)
                    else:
                        M_chi_r2[i].remove(tx1)    

    #the state before second iota, note the second round xor 0x2c
    for z in range(Zsize):
        con = 0x2c
        if ((con >> z) & 1) == 1:
            if ['1'] not in M_chi_r2[12*z]:
                M_chi_r2[12*z].append(['1'])
            else:
                M_chi_r2[12*z].remove(['1'])

    del M_chi_r3
    del M_roest_r2
    # state after the chi in the first round
    #only focus on planet[2] and planet[3]
    M_d = [[0 for i in range(vars_ln + 1)] for j in range(12 * Zsize)]
    for z in range(Zsize):
        for cell in range(6,12):
            i = cell + 12 * z
            for col in range(blocksize):
                if inv_Ma[i][col] == 1:
                    for term in M_chi_r2[col]:
                        if term == ['1']:
                            M_d[i][vars_ln] ^= 1
                        else:
                            p = Total_vars.index(term)
                            M_d[i][p] ^= 1
    del M_chi_r2
    #x_i = c_i + k_i
    A = [[0 for i in range(vars_ln + 1)] for j in range((n_out-3) * 3 * Zsize)]
    
    for m in range(n_out):
        #print('m = ' + str(m))
        C = [0 for i in range(blocksize)]
        for i in range(blocksize):
            x = int((i%12)/3)
            y = (i%12)%3
            z = int(i/12)
            C[i] = CC[m][x][y][z]  
        for z in range(Zsize):
            for cell in range(9,12):
                i = cell + 12 * z
                if m < (n_out-3):
                    flag = 3 * Zsize * m + 3 * z + (cell - 9)
                    #constant
                    A[flag][vars_ln] ^= M_d[i][vars_ln]                    
                    for k in range(vars_ln):
                        if M_d[i][k] == 1:
                            ln = len(Total_vars[k])
                            if ln == 1:
                                A[flag][k] ^= 1
                                A[flag][vars_ln] ^= C[Total_vars[k][0]]
                            
                            if ln == 2:
                                A[flag][k] ^= 1
                                p0 = Total_vars.index([Total_vars[k][0]])
                                p1 = Total_vars.index([Total_vars[k][1]])
                                A[flag][p0] ^= C[Total_vars[k][1]]
                                A[flag][p1] ^= C[Total_vars[k][0]]
                                A[flag][vars_ln] ^= (C[Total_vars[k][1]] * C[Total_vars[k][0]])
                                
                            if ln == 3:
                                #degree with 1
                                for ll in range(ln):
                                    p = Total_vars.index([Total_vars[k][ll]])
                                    A[flag][p] ^= (C[Total_vars[k][(ll + 1)%ln]] * C[Total_vars[k][(ll + 2)%ln]])
                                
                                #degree with 2
                                p = Total_vars.index([Total_vars[k][0], Total_vars[k][1]])
                                A[flag][p] ^= C[Total_vars[k][2]]
                                p = Total_vars.index([Total_vars[k][0], Total_vars[k][2]])
                                A[flag][p] ^= C[Total_vars[k][1]]
                                p = Total_vars.index([Total_vars[k][1], Total_vars[k][2]])
                                A[flag][p] ^= C[Total_vars[k][0]]
                                
                                #degree with 3
                                A[flag][k] ^= 1
                                
                                A[flag][vars_ln] ^= (C[Total_vars[k][0]] * C[Total_vars[k][1]] * C[Total_vars[k][2]])
                                
                            if ln == 4:
                                #degree with 1
                                for ll in range(ln):
                                    p = Total_vars.index([Total_vars[k][ll]])
                                    A[flag][p] ^= (C[Total_vars[k][(ll + 1)%ln]] * C[Total_vars[k][(ll + 2)%ln]] * C[Total_vars[k][(ll + 3)%ln]])
                                #degree with 2
                                for L1 in range(ln-1):
                                    for L2 in range(L1+1, ln):
                                        p = Total_vars.index([Total_vars[k][L1], Total_vars[k][L2]])
                                        tem = 1
                                        for re in range(ln):
                                            if (re != L1) and (re != L2):
                                                tem *= C[Total_vars[k][re]]
                                        A[flag][p] ^= tem
                                #degree with 3
                                for L1 in range(ln):
                                    list1 = []
                                    for L2 in range(ln):
                                        list1.append(Total_vars[k][L2])
                                    del list1[L1]
                                    p = Total_vars.index(list1)
                                    A[flag][p] ^= C[Total_vars[k][L1]]
                                #degree with 4
                                A[flag][k] ^= 1
                                
                                A[flag][vars_ln] ^= (C[Total_vars[k][0]] * C[Total_vars[k][1]] * C[Total_vars[k][2]]* C[Total_vars[k][3]])
                                
                                
                if m > 3:
                    i -= 3
                    flag = 3 * Zsize * (m-3) + 3 * z + (cell - 9)
                    #constant
                    A[flag][vars_ln] ^= M_d[i][vars_ln]                    
                    for k in range(vars_ln):
                        if M_d[i][k] == 1:
                            ln = len(Total_vars[k])
                            if ln == 1:
                                A[flag][k] ^= 1
                                A[flag][vars_ln] ^= C[Total_vars[k][0]]
                            
                            if ln == 2:
                                A[flag][k] ^= 1
                                p0 = Total_vars.index([Total_vars[k][0]])
                                p1 = Total_vars.index([Total_vars[k][1]])
                                A[flag][p0] ^= C[Total_vars[k][1]]
                                A[flag][p1] ^= C[Total_vars[k][0]]
                                
                                A[flag][vars_ln] ^= (C[Total_vars[k][1]] * C[Total_vars[k][0]])
                                
                            if ln == 3:
                                #degree with 1
                                for ll in range(ln):
                                    p = Total_vars.index([Total_vars[k][ll]])
                                    A[flag][p] ^= (C[Total_vars[k][(ll + 1)%ln]] * C[Total_vars[k][(ll + 2)%ln]])
                                
                                #degree with 2
                                p = Total_vars.index([Total_vars[k][0], Total_vars[k][1]])
                                A[flag][p] ^= C[Total_vars[k][2]]
                                p = Total_vars.index([Total_vars[k][0], Total_vars[k][2]])
                                A[flag][p] ^= C[Total_vars[k][1]]
                                p = Total_vars.index([Total_vars[k][1], Total_vars[k][2]])
                                A[flag][p] ^= C[Total_vars[k][0]]
                                
                                #degree with 3
                                A[flag][k] ^= 1
                                
                                A[flag][vars_ln] ^= (C[Total_vars[k][0]] * C[Total_vars[k][1]] * C[Total_vars[k][2]])
                                
                            if ln == 4:
                                #degree with 1
                                for ll in range(ln):
                                    p = Total_vars.index([Total_vars[k][ll]])
                                    A[flag][p] ^= (C[Total_vars[k][(ll + 1)%ln]] * C[Total_vars[k][(ll + 2)%ln]] * C[Total_vars[k][(ll + 3)%ln]])
                                #degree with 2
                                for L1 in range(ln-1):
                                    for L2 in range(L1+1, ln):
                                        p = Total_vars.index([Total_vars[k][L1], Total_vars[k][L2]])
                                        tem = 1
                                        for re in range(ln):
                                            if (re != L1) and (re != L2):
                                                tem *= C[Total_vars[k][re]]
                                        A[flag][p] ^= tem
                                #degree with 3
                                for L1 in range(ln):
                                    list1 = []
                                    for L2 in range(ln):
                                        list1.append(Total_vars[k][L2])
                                    del list1[L1]
                                    p = Total_vars.index(list1)
                                    A[flag][p] ^= C[Total_vars[k][L1]]
                                #degree with 4
                                A[flag][k] ^= 1
                                
                                A[flag][vars_ln] ^= (C[Total_vars[k][0]] * C[Total_vars[k][1]] * C[Total_vars[k][2]]* C[Total_vars[k][3]])                

    print('total vars = ' + str(vars_ln))
    flag = 0
    for i in range(vars_ln):
        s = 0
        for j in range((n_out-3) * 3 * Zsize):
            s += A[j][flag]
        if s == 0:
            #print(Total_vars[flag])
            del Total_vars[flag]
            for k in range ((n_out-3) * 3 * Zsize):
                del A[k][flag]
        else:
            flag += 1

    Candidate_key = solve_equations(A)
    
    if Candidate_key == K:
        print('Candidate Key:')
        print(Candidate_key)
        return 'success'
    else:
        return 'failed'
   

#3x+y+12z, write down the matrix of rho_{west}
Mp = [[0 for y in range(blocksize)] for x in range(blocksize)]
for i in range(blocksize):
    x = int((i%12)/3)
    y = (i%12)%3
    z = int(i/12)
    if y == 0:
        Mp[i][i] = 1

    if y == 1:
        Mp[i][3*((x-1)%4) + y + 12*z] = 1

    if y == 2:
        Mp[i][(i-var_rhowest*12)%blocksize] = 1
#3x+y+12z, write down the matrix of theta
Mt = [[0 for x in range(blocksize)] for y in range(blocksize)]
for x in range(Xsize):
    for y in range(Ysize):
        for z in range(Zsize):
            Mt[3*x+y+12*z][3*((x-1)%4) + 0 + 12*((z-var_theta1)%Zsize)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 1 + 12*((z-var_theta1)%Zsize)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 2 + 12*((z-var_theta1)%Zsize)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 0 + 12*((z-var_theta2)%Zsize)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 1 + 12*((z-var_theta2)%Zsize)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 2 + 12*((z-var_theta2)%Zsize)] = 1
            Mt[3*x+y+12*z][3*x+y+12*z] = 1
#3x+y+12z, write down the matrix of rho_{east}
Me = [[0 for y in range(blocksize)] for x in range(blocksize)]
for i in range(blocksize):
    x = int((i%12)/3)
    y = (i%12)%3
    z = int(i/12)
    if y == 0:
        Me[i][i] = 1

    if y == 1:
        Me[i][(i-1*12)%blocksize] = 1

    if y == 2:
        Me[i][3*((x-2)%4) + y + ((z-var_rhoeast)%Zsize)*12] = 1


Mb = multiply_binary_matrix(Mt,Me)
Ma = multiply_binary_matrix(Mp,Mb)

inv_Ma = inverse_binary_matrix(Ma)

result = three_round_attack(inv_Ma)
print(result)
            
            
            

    
