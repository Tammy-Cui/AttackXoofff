# write down this file on 25 April 2020
# this file is to verify the result of attack on 1-round xoofff



#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random



#define parameters
blocksize = 384
Xsize = 4
Ysize = 3
Zsize = 32
rd_xoodoo = 1 # total rounds of xoodoo used in xoofff
n_out = 15 # the number of output blocks in xoofff



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
        #A = rho_east_function(A)   # at the last round, we don't care this operation
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                A[x][y][z] =  A[x][y][z] ^ K[x][y][z]
    return A

        

def one_round_attack():
    
    S = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    C = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    K = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    # random choose S
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                S[x][y][z] = random.randint(0,1)
    
    for x in range(Xsize):
        for y in range(Ysize):
            for z in range(Zsize):
                K[x][y][z] = random.randint(0,1)
    #print(K)            
    #date stream is SS, output stream of xoofff is CC
    SS = [S for i in range(n_out)]
    CC = [C for i in range(n_out)]
    
    
    for n in range(n_out):
        CC[n] = Xoodoo_function(SS[n], K)
        if n < (n_out -1):
            SS[n+1] = roll_e_NLFSR_function(SS[n])

    
    k_guess = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    Candidate_K = [[[0 for z in range(Zsize)] for y in range(Ysize)] for x in range(Xsize)]
    for z in range(Zsize):
        candidate_k_num = 0
        #step 1
        #flag1 S'^{i}[2][2] = S'^{i+3}[1][2], and flag2 S'^{i}[2][0] = S'^{i+3}[1][0] after theta and rho_west        
        for k in range(1<<6):

            #guess 6 bit key k[1,y,z] and k[2,y,z]
            for y in range(Ysize):
                k_guess[1][y][z] = (k>>(2*y)) & 1
                k_guess[2][y][z] = (k>>(2*y+1)) & 1

            # n_out-3 conditions to filter the wrong 6-bit key
            test = 1
            for g in range(n_out-3):           
                flag1 = (CC[g][2][2][z] ^ k_guess[2][2][z]) ^ (CC[g][2][1][z] ^ k_guess[2][1][z]) ^ ((CC[g][2][0][z] ^ k_guess[2][0][z]) & (CC[g][2][1][z] ^ k_guess[2][1][z]))\
                        ^ (CC[g+3][1][2][z] ^ k_guess[1][2][z]) ^ (CC[g+3][1][1][z] ^ k_guess[1][1][z]) ^ ((CC[g+3][1][0][z] ^ k_guess[1][0][z]) & (CC[g+3][1][1][z] ^ k_guess[1][1][z]))
            
                if flag1 == 1:                    
                    test = 0
                    break
                
                flag2 = (CC[g][2][0][z] ^ k_guess[2][0][z]) ^ (CC[g][2][2][z] ^ k_guess[2][2][z]) ^ ((CC[g][2][1][z] ^ k_guess[2][1][z]) & (CC[g][2][2][z] ^ k_guess[2][2][z]))\
                        ^ (CC[g+3][1][0][z] ^ k_guess[1][0][z]) ^ (CC[g+3][1][2][z] ^ k_guess[1][2][z]) ^ ((CC[g+3][1][1][z] ^ k_guess[1][1][z]) & (CC[g+3][1][2][z] ^ k_guess[1][2][z]))
            
                if flag2 == 1:
                    test = 0
                    break
            
            if test == 1:
                
                #step 2
                #calculate k[3][y] form k[2][y], k[3][y] = k[2][y] ^ C^0[2][y] ^ C^3[3][y]                     
                for y in range(3):
                    k_guess[3][y][z] = k_guess[2][y][z] ^ CC[0][3][y][z] ^ CC[3][2][y][z]
                
                #step 3
                #flag3 S'^{i}[0][1] = S'^{i+3}[3][1] after theta and rho_west
                #guess 3-bit k[0,y,z]
                for k0 in range(1<<3):
                    for y in range(Ysize):
                        k_guess[0][y][z] = (k0>>y) & 1
                    tt = 1
                    for gg in range(n_out-3):
                        flag3 = (CC[gg][0][1][z] ^ k_guess[0][1][z]) ^ (CC[gg][0][0][z] ^ k_guess[0][0][z]) ^ ((CC[gg][0][2][z] ^ k_guess[0][2][z]) & (CC[gg][0][0][z] ^ k_guess[0][0][z]))\
                                ^(CC[gg+3][3][1][z] ^ k_guess[3][1][z]) ^ (CC[gg+3][3][0][z] ^ k_guess[3][0][z]) ^ ((CC[gg+3][3][2][z] ^ k_guess[3][2][z]) & (CC[gg+3][3][0][z] ^ k_guess[3][0][z]))
                        if flag3 == 1:
                            tt = 0
                            break
                    if tt == 1:
                        candidate_k_num += 1
                        for x in range(Xsize):
                            for y in range(Ysize):
                                Candidate_K[x][y][z] = k_guess[x][y][z]
                        

          
        if candidate_k_num > 1:
            return 'failed'
        
    #print(Candidate_K)
    return 'success'
    

### main function
succ = 0
T = 1000 #number of experiments
for t in range(T):
    flag = one_round_attack()
    if flag == 'success':
        succ += 1
print(str(succ/1000))


            
#A = Xoofff_function(A,K)
#print(A)

            
            
            

    
