import numpy as np
import random
from scipy.stats import expon

f = open('answer.txt', 'r+')
np.set_printoptions(suppress=True)

pn1 = 0.354
pn2 = 0.315
pn3 = 1 - pn1 - pn2
pm0 = 0.525
pm1 = 0.299
pm2 = 1 - pm0 - pm1

N = 0
M = 0

L = list(range(1, 101))

Ttime = []
Ttype = []
Tl1 = []
Tl2 = []
condition = []
Tremained = []
Tdin = []
TdG = []

jN = np.zeros(300)
jG = np.zeros(300)
jB = np.zeros(300)
jLt = np.zeros(300)
jD = np.zeros(300)
jS1 = np.zeros(300)
jS2 = np.zeros(300)

Sn1 = 0
Sn2 = 0
Sn3 = 0

Sm0 = 0
Sm1 = 0
Sm2 = 0

Sn1d = 0
Sn2d = 0
Sn3d = 0
Sm0d = 0
Sm1d = 0
Sm2d = 0


n5 = [[0 for j in range(15)] for j in range(50)]
v5 = [[0 for j in range(15)] for j in range(50)]
T5 = [[0 for j in range(15)] for j in range(50)]
D5 = [[0 for j in range(15)] for j in range(50)]

Ld = []
Ldi = []

arrive_time = expon.rvs(scale=1 / 1)

i = 1

Ttime.append(0)
Ttype.append("Sn1")
Tl1.append(expon.rvs(scale=1 / 0.3))
N += 1
Tl2.append(-1)
condition.append((1, 0))
Ld.append(Tl1[0])
Ldi.append("N")
Tremained.append(Tl1[0])
Tdin.append(i)
TdG.append("N")

jN[0] = i
jG[0] = 1
jB[0] = 0
jLt[0] = Tl1[0]
jD[0] = Tl1[0]

Sn1 += 1
n5[1][0] += 1
T5[1][0] += Tl1[0]
# T1 += service_time[i]

while len(Ttime) != 100:
    if Ldi[np.argmin(Ld)] == "N":

        Ttime.append(Ttime[-1] + Ld[np.argmin(Ld)])
        с = random.random()
        if с < pn1:
            jS1[np.argmin(Ld)] = i + 1
            jS2[np.argmin(Ld)] = -1

            Ttype.append("Sn1")
            Tl1.append(expon.rvs(scale=1 / (0.005 * N + 0.01 * M + 0.3)))
            Tl2.append(-1)
            condition.append(condition[-1])
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))
            Ld = np.append(Ld, Tl1[-1])
            Ldi.append("N")
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld)+1)
            TdG.append(Ldi[np.argmin(Ld)])



            jN[i] = i
            jG[i] = 1
            jB[i] = Ttime[-1]
            jLt[i] = Ld[i]
            jD[i] = Ttime[-1] + Ld[i]
            i += 1
            Sn1 += 1
            n5[N][M] = n5[N][M]+ 1
            T5[N][M] += Tremained[-1]
        elif с < pn2 + pn1:
            jS1[np.argmin(Ld)] = i + 1
            jS2[np.argmin(Ld)] = i + 2

            Ttype.append("Sn2")
            a = expon.rvs(scale=1 / (0.005 * N + 0.01 * M + 0.3))
            b = expon.rvs(scale=1 / (0.005 * N + 0.01 * M + 0.3))
            if b > a:
                Tl1.append(a)
                Tl2.append(b)
            else:
                Tl1.append(b)
                Tl2.append(a)
            condition.append((condition[-1][0] + 1, condition[-1][1]))
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))
            Ld = np.append(Ld, Tl1[-1])
            Ldi.append("N")
            Ld = np.append(Ld, Tl2[-1])
            Ldi.append("N")
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld)+1)
            TdG.append(Ldi[np.argmin(Ld)])



            jN[i] = i + 1
            jG[i] = 1
            jB[i] = Ttime[-1]
            jLt[i] = Ld[i]
            jD[i] = Ttime[-1] + Ld[i]
            i += 1

            jN[i] = i + 1
            jG[i] = 2
            jB[i] = Ttime[-1]
            jLt[i] = Ld[i]
            jD[i] = Ttime[-1] + Ld[i]
            i += 1

            Sn2 += 1
            N += 1
            n5[N][M] += 1
            T5[N][M] += Tremained[-1]
        else:
            Ttype.append("Sn3")
            Sn3+=1

            jS1[np.argmin(Ld)] = i + 1
            jS2[np.argmin(Ld)] = i + 2

            a = expon.rvs(scale=1 / (0.005 * N + 0.01 * M + 0.3))
            b = expon.rvs(scale=1 / (0.05 * N + 0.03 * M))
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))


            if b > a:
                Tl1.append(a)
                Tl2.append(b)
                Ld = np.append(Ld, Tl1[-1])
                Ldi.append("N")
                Ld = np.append(Ld, Tl2[-1])
                Ldi.append("M")

                jN[i] = i + 1
                jG[i] = 1
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

                jN[i] = i + 1
                jG[i] = 2
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1
            else:
                Tl1.append(b)
                Tl2.append(a)
                Ld = np.append(Ld, (Tl1[-1]))
                Ldi.append("M")
                Ld = np.append(Ld, (Tl2[-1]))
                Ldi.append("N")

                jN[i] = i + 1
                jG[i] = 2
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

                jN[i] = i + 1
                jG[i] = 1
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

            condition.append((condition[-1][0], condition[-1][1] + 1))
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld)+1)
            TdG.append(Ldi[np.argmin(Ld)])
            M += 1
            n5[N][M] += 1
            T5[N][M] += Tremained[-1]
    else:
        Ttime.append(Ttime[-1] + Ld[np.argmin(Ld)])
        с = random.random()
        if с < pn1:
            Ttype.append("Sm1")

            jS1[np.argmin(Ld)] = i + 1
            jS2[np.argmin(Ld)] = -1

            Tl1.append(expon.rvs(scale=1 / (0.05 * N + 0.03 * M)))
            Tl2.append(-1)
            condition.append(condition[-1])
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))
            Ld = np.append(Ld, (Tl1[-1]))
            Ldi.append("M")
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld) + 1)
            TdG.append(Ldi[np.argmin(Ld)])



            jN[i] = i + 1
            jG[i] = 2
            jB[i] = Ttime[-1]
            jLt[i] = Ld[i]
            jD[i] = Ttime[-1] + Ld[i]
            i += 1

            Sm1 += 1
            n5[N][M] += 1
            T5[N][M] += Tremained[-1]
        elif с < pn2 + pn1:
            Ttype.append("Sm0")

            jS1[np.argmin(Ld)] = -1
            jS2[np.argmin(Ld)] = -1

            Tl1.append(-1)
            Tl2.append(-1)
            condition.append((condition[-1][0], condition[-1][1] - 1))
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld)+1)
            TdG.append(Ldi[np.argmin(Ld)])


            Sm0 += 1
            M -= 1
            n5[N][M] += 1
            T5[N][M] += Tremained[-1]
        else:
            Ttype.append("Sm2")

            jS1[np.argmin(Ld)] = i + 1
            jS2[np.argmin(Ld)] = i + 2

            a = expon.rvs(scale=1 / (0.005 * N + 0.01 * M + 0.3))
            b = expon.rvs(scale=1 / (0.05 * N + 0.03 * M))
            Ld -= Ld[np.argmin(Ld)]
            Ld[np.argmin(Ld)] = (float("inf"))


            if b > a:
                Tl1.append(a)
                Tl2.append(b)
                Ld = np.append(Ld, Tl1[-1])
                Ldi.append("N")
                Ld = np.append(Ld, (Tl2[-1]))
                Ldi.append("M")

                jN[i] = i + 1
                jG[i] = 1
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

                jN[i] = i + 1
                jG[i] = 2
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1
            else:
                Tl1.append(b)
                Tl2.append(a)
                Ld = np.append(Ld, (Tl1[-1]))
                Ldi.append("M")
                Ld = np.append(Ld, (Tl2[-1]))
                Ldi.append("N")

                jN[i] = i + 1
                jG[i] = 2
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

                jN[i] = i + 1
                jG[i] = 1
                jB[i] = Ttime[-1]
                jLt[i] = Ld[i]
                jD[i] = Ttime[-1] + Ld[i]
                i += 1

            condition.append((condition[-1][0]+1, condition[-1][1]))
            Tremained.append(Ld[np.argmin(Ld)])
            Tdin.append(np.argmin(Ld)+1)
            TdG.append(Ldi[np.argmin(Ld)])
            Sm2 += 1
            N += 1
            n5[N][M] += 1
            T5[N][M] += Tremained[-1]

Sn1d = Sn1 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)
Sn2d = Sn2 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)
Sn3d = Sn3 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)
Sm0d = Sm0 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)
Sm1d = Sm1 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)
Sm2d = Sm2 / (Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2)

T5[N][M] -= Tremained[-1]
v5 = np.array(n5) / 100
D5 = np.array(T5) / Ttime[-1]

Gn = 0
Gm = 0
for i in range(len(Ld)):
    if (Ldi[i] == "N"):
        Gn += 1
    else:
        Gm += 1

f.write(str(np.around(L, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttime, 5)))
f.write('\n')
f.write('\n')
f.write(str(Ttype))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Tl1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tl2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(condition, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tremained, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tdin, 5)))
f.write('\n')
f.write('\n')
f.write(str(TdG))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jN, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jG, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jB, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jLt, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jD, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS2, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn3, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm0, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn1d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn2d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sn3d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm0d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm1d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Sm2d, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around((Sn1 + Sn2 + Sn3 + Sm1 + Sm0 + Sm2), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around((Sn1d + Sn2d + Sn3d + Sm1d + Sm0d + Sm2d), 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Gn, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Gm, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(N, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(M, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(n5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(v5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(T5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(D5, 5)))
f.write('\n')
f.write('\n')
total = 0
for ele in range(0, 50):
    total += sum(n5[ele])
f.write(str(np.around(total, 5)))
f.write('\n')
f.write('\n')
total = 0
for ele in range(0, 50):
    total += sum(v5[ele])
f.write(str(np.around(total, 5)))
f.write('\n')
f.write('\n')
total = 0
for ele in range(0, 50):
    total += sum(T5[ele])
f.write(str(np.around(total, 5)))
f.write('\n')
f.write('\n')
total = 0
for ele in range(0, 50):
    total += sum(D5[ele])
f.write(str(np.around(total, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
