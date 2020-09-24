# write down this file on 25 April 2020
# this file is to verify the result of attacks on 2-round xoofff



#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random



#define parameters
blocksize = 384
Xsize = 4
Ysize = 3
Zsize = 32
rd_xoodoo = 2 # total rounds of xoodoo used in xoofff
n_out = 100 # the number of output blocks of xoofff
n_step1 = 18
n_step2 = 40
n_step3 = 73
n_step4 = 6

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
    #print('row = ' + str(row))
    #print('col = ' + str(col))
    if row == 1:
        for i in range(col):
            if B[row-1][i] == 1:
                return 1
        return 0
    #col-1 * col-1 matrix
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
            Guess_K = [0 for i in range(col-1)]
            
            for i in range(col-1):
                if B[col-2-i][col-2-i] == 0:
                    #print('k'+str(col-2-i) +' is free variable')
                    return 'failed'
                if i == 0:
                    Guess_K[col - 2 - i] = B[col-2][col - 1]
                else:
                    Guess_K[col - 2 - i] = B[col-2-i][col - 1]
                    for j in range(i):
                        Guess_K[col - 2 - i] ^= B[col-2-i][col-2-j] * Guess_K[col - 2 - j]
            #print(Guess_K)
            
            return Guess_K
def inverse_binary_matrix(A):
    row = len(A)
    col = len(A[0])
    if row != col:
        return "failed"
    #rank = rank_binary_matrix(A)
    #if rank != row:
     #   return "error"
    
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
            
    # E = P <<< (1,5) + P <<< (1,14)
    E = [[0 for z in range(Zsize)] for x in range(Xsize)]
    for x in range(Xsize):
        for z in range(Zsize):
            E[x][z] = P[(x-1) % Xsize][(z-5) % Zsize] ^ P[(x-1) % Xsize][(z-14) % Zsize]

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
            #A_2 = A_2 <<< (0,11)
            C[x][z] = A[x][2][(z-11) % Zsize]

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
            #A_2 = A_2 <<< (2,8)
            C[x][z] = A[(x-2) % Xsize][2][(z-8) % Zsize]

    for x in range(Xsize):
        for z in range(Zsize):
            A[x][1][z] = B[x][z]
            A[x][2][z] = C[x][z]
    return A

# NLFSR, roll_e
def roll_e_NLFSR_function(A):
    
    #A[0,0] = A[0,1]* A[0,2] + A[0,0]<<<5 + A[0,1]<<<13 + 0x00000007
    C = [0 for z in range(Zsize)]
    for z in range(Zsize):
        C[z] = (A[0][1][z] & A[0][2][z]) ^ A[0][0][(z-5) % Zsize] ^ A[0][1][(z-13) % Zsize] ^ ((0x00000007 >> z) & 1)

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


def two_round_attack():
    
    S = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    C = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    K = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    # random choose input state S of expansion layer and masterkey K
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                S[x][y][z] = random.randint(0,1)

    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                K[x][y][z] = random.randint(0,1)
    #print("True Key:")
    #print(K)
                    
    #date stream is SS, output stream of xoofff is CC
    SS = [S for i in range(n_out)]
    CC = [C for i in range(n_out)]
    
    #o = open("test.txt", 'w')
    for n in range(n_out):
        CC[n] = Xoodoo_function(SS[n], K)
        if n < (n_out -1):
            SS[n+1] = roll_e_NLFSR_function(SS[n])
    
    #step 1
    k_guess = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]

     
    for z in range(Zsize):
        for k in range(1<<9):
            #guess 9 bit key k[0,y,z] and k[2,y,z], k[3,y,z]
            for y in range(Ysize):
                k_guess[0][y][z] = (k>>(3*y)) & 1
                k_guess[2][y][z] = (k>>(3*y+1)) & 1
                k_guess[3][y][z] = (k>>(3*y+2)) & 1

            # n_out-3 conditions to filter the wrong 6-bit key
            test = 1
            for g in range(n_step1-3):
                flag1 = (CC[g][0][1][z] ^ k_guess[0][1][z]) ^ (CC[g][0][0][z] ^ k_guess[0][0][z]) ^ ((CC[g][0][2][z] ^ k_guess[0][2][z]) & (CC[g][0][0][z] ^ k_guess[0][0][z]))\
                        ^ (CC[g][3][0][z] ^ k_guess[3][0][z]) ^ (CC[g][3][2][z] ^ k_guess[3][2][z]) ^ ((CC[g][3][1][z] ^ k_guess[3][1][z]) & (CC[g][3][2][z] ^ k_guess[3][2][z]))\
                        ^ (CC[g+3][3][1][z] ^ k_guess[3][1][z]) ^ (CC[g+3][3][0][z] ^ k_guess[3][0][z]) ^ ((CC[g+3][3][2][z] ^ k_guess[3][2][z]) & (CC[g+3][3][0][z] ^ k_guess[3][0][z]))\
                        ^ (CC[g+3][2][0][z] ^ k_guess[2][0][z]) ^ (CC[g+3][2][2][z] ^ k_guess[2][2][z]) ^ ((CC[g+3][2][1][z] ^ k_guess[2][1][z]) & (CC[g+3][2][2][z] ^ k_guess[2][2][z]))
                    
                if flag1 == 1:
                    test = 0
                    break
                
            
            if test == 1:
                break



    #step 2: check 2^32 possible cases (k[0,1],k[2,0]) by parity on MMx = BB
    MM = [[0 for y in range(32+1)] for x in range(n_step2 - 3)]
    #BB = [0 for x in range(n_step2 - 3)]

    for i in range(n_step2 - 3):        
        for z in range(Zsize):
            MM[i][z] = CC[i+3][2][2][z] ^ k_guess[2][2][z] ^ 1
            MM[i][32] ^= (CC[i][3][1][z] ^ k_guess[3][1][z]) ^ ((CC[i][3][2][z] ^ k_guess[3][2][z] ^ 1) & (CC[i][3][0][z] ^ k_guess[3][0][z]))
            MM[i][32] ^= (CC[i+3][2][1][z] ^ k_guess[2][1][z]) ^ ((CC[i+3][2][2][z] ^ k_guess[2][2][z] ^ 1) & CC[i+3][2][0][z])



    var_x = solve_equations(MM)
    if var_x == 'failed':
        return 'failed'

    for z in range(Zsize):
        if var_x[z] != k_guess[2][0][z]:
            k_guess[2][0][z] ^= 1
            k_guess[0][1][z] ^= 1


    #until now, we have recover 288-bit key on sheet 0, 2 and 3.
  
    #step 3: recover 64-bit key k[1,0] and k[1,1]
    #Mx = B
    M = [[0 for y in range(64+1)] for x in range(n_step3 - 4)]
    #B = [0 for x in range((n_step3 - 3)*(n_step3 - 3))]

    for i in range(n_step3-4):
        j = i + 1
        for z in range(Zsize):
            M[i][z] = CC[i][1][1][z] ^ CC[j][1][1][z]
            M[i][32+z] = CC[i][1][0][z] ^ CC[j][1][0][z]
                    
            M[i][64] ^= CC[i][1][2][z] ^ CC[i][1][1][z] ^ (CC[i][1][0][z] & CC[i][1][1][z])
            M[i][64] ^= CC[j][1][2][z] ^ CC[j][1][1][z] ^ (CC[j][1][0][z] & CC[j][1][1][z])
            M[i][64] ^= (CC[i+3][0][2][z] ^k_guess[0][2][z]) ^ (CC[i+3][0][1][z] ^k_guess[0][1][z]) ^ ((CC[i+3][0][0][z] ^k_guess[0][0][z]) & (CC[i+3][0][1][z] ^k_guess[0][1][z]))
            M[i][64] ^= (CC[j+3][0][2][z] ^k_guess[0][2][z]) ^ (CC[j+3][0][1][z] ^k_guess[0][1][z]) ^ ((CC[j+3][0][0][z] ^k_guess[0][0][z]) & (CC[j+3][0][1][z] ^k_guess[0][1][z]))

    var_x = solve_equations(M)
    if var_x == 'failed':
        return 'failed'
                    
    for z in range(Zsize):
        k_guess[1][0][z] = var_x[z]
        k_guess[1][1][z] = var_x[Zsize + z]
   

    #step4: recover last 32-bit key k[1,2]     

    #Mx = B
    M = [[0 for y in range(32+1)] for x in range((n_step4-3) * Zsize)]
    #B = [0 for y in range(32)]
    n = 0
    for dt in range(n_step4-3):
        for eq in range(Zsize):
            for i in range(blocksize):
                x = int((i%12)/3)
                y = (i%12)%3
                z = int(i/12)

                M[n][32] ^= inv_Ma[10+12*eq][i] & ((CC[dt][x][y][z] ^ k_guess[x][y][z]) ^ ((CC[dt][x][(y+1)%3][z] ^ k_guess[x][(y+1)%3][z] ^ 1) & (CC[dt][x][(y+2)%3][z] ^ k_guess[x][(y+2)%3][z])))
                M[n][32] ^= inv_Ma[7+12*eq][i] &  ((CC[dt+3][x][y][z] ^ k_guess[x][y][z]) ^ ((CC[dt+3][x][(y+1)%3][z] ^ k_guess[x][(y+1)%3][z] ^ 1) & (CC[dt+3][x][(y+2)%3][z] ^ k_guess[x][(y+2)%3][z])))
                if i == 24 or i == 36 or i == 60:
                    M[n][32] ^= inv_Ma[10+12*eq][i] & 1
                    M[n][32] ^= inv_Ma[7+12*eq][i] & 1
                    
            for var in range(32):
                M[n][var] ^= inv_Ma[10+12*eq][3+12 * var] & (CC[dt][1][1][var] ^ k_guess[1][1][var] ^ 1)
                M[n][var] ^= inv_Ma[10+12*eq][4+12 * var] & (CC[dt][1][0][var] ^ k_guess[1][0][var])
                M[n][var] ^= inv_Ma[10+12*eq][5+12 * var] & 1

                M[n][var] ^= inv_Ma[7+12*eq][3+12 * var] & (CC[dt+3][1][1][var] ^ k_guess[1][1][var] ^ 1)
                M[n][var] ^= inv_Ma[7+12*eq][4+12 * var] & (CC[dt+3][1][0][var] ^ k_guess[1][0][var])
                M[n][var] ^= inv_Ma[7+12*eq][5+12 * var] & 1            


            n += 1

    var_x = solve_equations(M)
    if var_x == 'failed':
        return 'failed'
    
    for z in range(Zsize):
        k_guess[1][2][z] = var_x[z]   

    if k_guess == K:
        #print("Candidate Key:")
        #print(k_guess)
        return "success"
    else:
        return "failed"



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
        Mp[i][(i-11*12)%384] = 1

#3x+y+12z, write down the matrix of theta
Mt = [[0 for x in range(blocksize)] for y in range(blocksize)]
for x in range(Xsize):
    for y in range(Ysize):
        for z in range(Zsize):
            Mt[3*x+y+12*z][3*((x-1)%4) + 0 + 12*((z-5)%32)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 1 + 12*((z-5)%32)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 2 + 12*((z-5)%32)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 0 + 12*((z-14)%32)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 1 + 12*((z-14)%32)] = 1
            Mt[3*x+y+12*z][3*((x-1)%4) + 2 + 12*((z-14)%32)] = 1
            Mt[3*x+y+12*z][3*x+y+12*z] = 1

Ma = multiply_binary_matrix(Mp,Mt)
inv_Ma = inverse_binary_matrix(Ma)





T = 100   #number of experiments
s = 0
for i in range(T):
    r = two_round_attack()
    if r == "success":
        s += 1
print (str(s/T))


            
            
            

    
