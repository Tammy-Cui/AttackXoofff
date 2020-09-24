# AttackXoofff

There are four files of practical attacks:

*PracticalAttack-OneRound-Xoofff.py
*PracticalAttack-TwoRound-Xoofff.py
*PracticalAttack-TwoRound-Xoofff(linearization).py
*PracticalAttack-ThreeRound-Xoofff.py

All practical attacks on Xoofff are aimed to recover the key under any random secret true key.

File 1: PracticalAttack-OneRound-Xoofff.py

Goal: to find out the success rate of key-recovery attack on 1-round Xoofff under different numbers of output blocks.

Modifiable variables: 
1, n_out: the number of output blocks used in the attack;
2, T: number of experiments

If you want to check the true secret key and candidate key, you can reopen the remarks:
#print(K)
#print(Candidate_K)
in function one_round_attack()


File 2: PracticalAttack-TwoRound-Xoofff.py

Goal: to find out the success rate of key-recovery attack on 2-round Xoofff under different numbers of output blocks.

Modifiable variables: 
1, n_out: the number of output blocks used in the attack;
2, n_step1, n_step2, n_step3, n_step4: the number of output blocks used in each step of attack;
3, T: number of experiments

File 3: PracticalAttack-TwoRound-Xoofff(linearization).py

Goal:  to recover the secret key with linearization technique on 2-round Xoofff (toy-version cipher)

There are three versions of cipher according to the length of lane, i.e. Zsize 4, 8, 32. We implement the attack with a toy-version cipher Zsize = 8 because of the time complexity.

Modifiable variables: 
1, n_out: the number of output blocks used in the attack;
2, Zsize: the length of lane, related to the version of cipher we use

File 4: PracticalAttack-ThreeRound-Xoofff.py

Goal:  to recover the secret key with linearization technique on 3-round Xoofff (toy-version cipher)

There are three versions of cipher according to the length of lane, i.e. Zsize 4, 8, 32. We implement the attack with a toy-version cipher Zsize =4 or 8 because of the time complexity.

Modifiable variables: 
1, n_out: the number of output blocks used in the attack;
2, Zsize: the length of lane, related to the version of cipher we use
