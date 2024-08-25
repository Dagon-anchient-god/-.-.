import math
import numpy as np
from scipy.stats import expon
from decimal import Decimal

f = open('answer.txt', 'r+')

np.set_printoptions(suppress=True)
u = 1.214
lambd = 5.865
Tt = 0.17
n = 5
m = 14
# Время между заявками

L = list(range(1, 101))
k = [float('inf'), float('inf'), float('inf'), float('inf'), float('inf')]
kj = [-1, -1, -1, -1, -1]
service_time = expon.rvs(scale=1 / u, size=100)
service_time = [0.2089425, 0.24042362, 0.17402998, 0.28750913, 0.57853876, 0.3281805, 0.64220276, 0.05176051, 0.51798679,
                0.12153939, 0.08446545, 0.51658193, 0.0755792, 0.96452253, 1.24035621, 0.1616758, 0.44373353, 0.02876248,
                0.10613567, 1.33902265, 1.01999454, 1.1779021, 0.02906335, 0.99077037, 0.68277089, 2.99662216, 1.16265632,
                0.24735046, 0.03655438, 0.63890711, 0.04985124, 3.34427941, 3.2353478, 1.65940145, 0.37414538, 0.43366022,
                0.10008147, 0.68081293, 0.11424763, 0.55604128, 0.30231069, 0.71615587, 1.31168242, 0.0363553, 0.91358345,
                0.43515835, 0.35815254, 0.07650709, 0.85657331, 0.03259491, 0.41606914, 2.12135799, 1.19049445, 0.11147616,
                0.64662056, 0.14681455, 0.65388675, 1.46099198, 2.28533155, 0.63652264, 0.03639166, 1.06086346, 0.20777427,
                0.64148374, 0.42320609, 0.55148203, 0.65518179, 0.46765486, 0.70194748, 0.1279088, 1.17681095, 0.3871591,
                1.23335323, 0.76315801, 0.90589213, 2.12009595, 0.25639081, 0.68965077, 1.14325731, 0.64423279, 1.22593882,
                0.48603955, 0.36128709, 0.71295673, 1.04207984, 0.7153235, 0.33775967, 0.66214228, 0.32824737, 2.90752264,
                1.37187689, 1.41407363, 0.82795116, 0.35135893, 0.03515003, 0.23388852, 0.55634401, 0.04169105, 0.07508141,
                0.24626046]
Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []
numk = []

jN = np.zeros(100)
jP = np.zeros(100)
jQ = np.zeros(100)
jQt = np.zeros(100)
jS = np.zeros(100)
jD = np.zeros(100)
jF = np.zeros(100)
jk = np.zeros(100)

R = np.zeros(100)
V = np.zeros(100)

k41 = np.zeros(n)
k42 = np.zeros(n)
k43 = np.zeros(n)
k44 = np.zeros(n)

i = 0
it = 0

S_con = 0
trimen = 0

R[0] += 1
V[0] += Tt

J5 = 0
JF5 = 0
Z5 = 0
Tq = 0
Tl = 0

Ttime.append(Tt)

jN[it] = i + 1
jP[it] = Tt
jQ[it] = 0
jQt[it] = 0
jS[it] = Tt
jD[it] = service_time[0]
jF[it] = Tt + service_time[0]
jk[it] = it + 1

Ttype.append(1)
condition.append(1)

# T1 += service_time[i]

Tremained.append(service_time[0])
Tnew.append(Tt)
numj.append(it + 1)
numk.append(it + 1)
trimen = Tt - service_time[0]
S_con = 1
k[0] = service_time[it]
kj[0] = it+1

i += 1
it += 1
J5 += 1


while len(Ttime) != 100:
    if S_con == 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + Tnew[-1]
        jQ[i] = 0
        jQt[i] = 0
        jS[i] = Ttime[-1] + Tnew[-1]
        jk[it] = S_con + 1

        V[S_con] += Tnew[-1]
        Ttime.append(Ttime[-1] + Tnew[-1])
        Ttype.append(1)
        condition.append(S_con + 1)

        # T1 += service_time[i]

        k[S_con] = service_time[i]
        Tremained.append(service_time[i])
        Tnew.append(Tt)
        numj.append(i + 1)
        numk.append(S_con + 1)
        trimen = Tt - service_time[i]
        S_con += 1
        R[S_con] += 1
        k[0] = service_time[i]
        kj[0] = i+1
        k41[0] += 1
        k42[0] += service_time[i]
        # V[S_con] += min(Ts, arrive_time[it + 1])
        it += 1
        i += 1
        J5 += 1
        Z5 += S_con
    elif trimen > 0:
        jS[kj[np.argmin(k)]-1] = Ttime[-1] + k[np.argmin(k)] - service_time[kj[np.argmin(k)]-1]
        jD[kj[np.argmin(k)]-1] = service_time[kj[np.argmin(k)]-1]
        jF[kj[np.argmin(k)]-1] = Ttime[-1] + Tremained[-1]
        jk[kj[np.argmin(k)]-1] = np.argmin(k) + 1

        V[S_con] += k[np.argmin(k)]
        numj.append(kj[np.argmin(k)])
        Ttime.append(Ttime[-1] + k[np.argmin(k)])
        k = np.subtract(k, k[np.argmin(k)])
        Ttype.append(2)
        condition.append(S_con - 1)
        if S_con <= n:
            kj[np.argmin(k)] = -1
            k[np.argmin(k)] = float('inf')
        else:
            k[np.argmin(k)] = float('inf')
            iop = np.argmax(k)
            jS[j] = Ttime[-1]
            jD[j] = service_time[j]
            jF[j] = Ttime[-1] + service_time[j]
            jk[j] = np.argmax(k) + 1
            k41[iop] += 1
            k42[iop] += service_time[kj[iop] - 1]
            k[iop] = service_time[j]
            kj[iop] = j+1
            j += 1
            if j>i:
                i+=1
        if S_con - 1 == 0:
            Tremained.append(-1)
            Tnew.append(trimen)
            Tl += trimen
        else:
            Tremained.append(min(k))
            Tnew.append(trimen)

        trimen = trimen - min(k)
        numk.append(np.argmin(k) + 1)
        S_con -= 1
        R[S_con] += 1
        JF5 += 1
        Z5 += S_con
    elif trimen < 0:
        tp = trimen + min(k)
        k = np.subtract(k, tp)
        jN[i] = i + 1
        jP[i] = Ttime[-1] + tp
        jQ[i] = S_con
        if S_con < n:
            jQt[i] = 0
        else:
            for d in range(S_con):
                jQt[i] += service_time[d + it + 1]
            jQt[i] -= trimen

        if max(k) == float('inf'):
            numk.append(np.argmax(k) + 1)
            jS[it] = (i + 1) * Tt + jQt[i]
            jD[it] = service_time[i]
            jF[it] = (i + 1) * Tt + jQt[i] + service_time[i]
            jk[i] = np.argmax(k) + 1
            k41[np.argmax(k)] += 1
            k42[np.argmax(k)] += service_time[kj[np.argmax(k)] - 1]
            kj[np.argmax(k)] = i+1
            k[np.argmax(k)] = service_time[i]
            j = i+1
        else:
            numk.append(-1)
        V[S_con] += tp
        Ttime.append((i + 1) * Tt)
        Ttype.append(1)
        condition.append(S_con + 1)
        Tremained.append(min(k))
        Tnew.append(Tt)
        i += 1
        S_con += 1
        R[S_con] += 1
        numj.append(i)
        trimen = Tt - min(k)
        J5 += 1
        Z5 += S_con

Rot = R / 100
Vot = V / Ttime[-1]
Z5 = Z5 / 100
wer = 0
for d3 in range(len(jF)):
    if jF[d3] > 0:
        wer += (jF[d3] - jS[d3])
Tq5 = sum(jQt) / JF5
Tm5 = wer / JF5

for d1 in range(len(jQ)):
    if jQ[d1] < 5:
        jQ[d1] = 0
    else:
        jQ[d1] -= 4

for d2 in range(n):
    if k[d2] != float('inf'):
        k42[d2] -= k[d2]
    k43[d2] = Ttime[-1] - k42[d2]
    k44[d2] = k43[d2] / Ttime[-1]

f.write(str(np.around(L, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttime, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttype, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(condition, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tremained, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tnew, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(numj, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(numk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jN, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jP, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQ, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQt, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jD, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jF, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(R, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(V, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Rot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Vot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(R), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(V), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Rot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Vot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(J5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(JF5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Z5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tq5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Tm5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(k41, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k42, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k43, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k44, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')

k = [float('inf'), float('inf'), float('inf'), float('inf'), float('inf')]
kj = [-1, -1, -1, -1, -1]
Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []
numk = []

jN = np.zeros(100)
jP = np.zeros(100)
jQ = np.zeros(100)
jQt = np.zeros(100)
jS = np.zeros(100)
jD = np.zeros(100)
jF = np.zeros(100)
jk = np.zeros(100)

R = np.zeros(100)
V = np.zeros(100)

k41 = np.zeros(n)
k42 = np.zeros(n)
k43 = np.zeros(n)
k44 = np.zeros(n)

T0 = Tt
T1 = 0

arrive_time = expon.rvs(scale=1 / lambd, size=105)
arrive_time = [0.13452634, 0.15519843, 0.0283299, 0.04746884, 0.05783817, 0.09628136, 0.16774346, 0.04122172, 0.26758845,
               0.02089449, 0.0062113, 0.1203606, 0.01673303, 0.49141359, 0.13555423, 0.24224803, 0.0384094, 0.01127594,
               0.20048834, 0.2155568, 0.33837971, 0.04805961, 0.08154375, 0.03341779, 0.01885529, 0.58151262, 0.18703896,
               0.03363032, 0.16354833, 1.05523772, 0.20401247, 0.0897726, 0.07839006, 0.00064297, 0.11498224, 0.0968871,
               0.18283561, 0.28195189, 0.09633881, 0.05745809, 0.29990406, 0.31262433, 0.52178662, 0.07495145, 0.32229107,
               0.69479259, 0.3485287, 0.00978547, 0.34455846, 0.14409691, 0.0015866, 0.24095841, 0.20724527, 0.14153678, 0.21186249,
               0.11534158, 0.27668605, 0.23990952, 0.12584108, 0.02072617, 0.18680627, 0.10086884, 0.23286964, 1.29150504,
               0.19749515, 0.14165763, 0.01949458, 0.07271728, 0.06531828, 0.16581318, 0.08882034, 0.26548143, 0.17186677,
               0.04648546, 0.06211913, 0.06518982, 0.12855165, 0.35920409, 0.0315275, 0.0637701, 0.23644445, 0.0669599,
               0.13345968, 0.09329895, 0.041594, 0.19948896, 0.20128085, 0.0413679, 0.05342318, 0.20120517, 0.3631961,
               0.38923137, 0.09698655, 0.04347043, 0.01412543, 0.38276278, 0.20002275, 0.07304276, 0.07767907, 0.20347489,
               0.0347942, 0.16591234, 0.12607382, 0.00059887, 0.18488854]
service_time = expon.rvs(scale=1 / u, size=105)
service_time = [0.14096427, 1.52599583, 1.36619083, 2.78187587, 0.01550708, 0.50005387, 0.31878136, 2.42388767, 0.47219829,
                0.34627476, 0.63164567, 0.31118265, 0.89021229, 0.57634907, 0.41271947, 0.34772789, 0.75406021, 0.40428931,
                0.15281013, 0.96501346, 0.01430193, 0.22747386, 0.42769759, 1.44223716, 0.1218457, 1.55408956, 0.25950125,
                1.14195904, 1.59255417, 0.2110721, 2.21962061, 0.03299088, 0.16090538, 0.85345717, 1.07927905, 1.24191927,
                0.69821569, 0.47136381, 0.85132492, 0.14742841, 0.57725097, 0.12689934, 0.24906341, 0.08551547, 0.86465192,
                0.08699023, 0.45010711, 1.72123805, 0.91003177, 0.86917588, 0.81969871, 1.00404901, 0.35450958, 0.13491099,
                0.61593932, 0.08734727, 1.32834542, 2.77345115, 2.84737523, 0.11407401, 2.18020131, 0.58461787, 1.20795147,
                0.17302712, 0.62418966, 0.29549693, 0.19365026, 0.12110518, 0.0660802, 0.73917226, 0.69735931, 0.37485889,
                0.5456252, 1.83985775, 2.70396212, 0.07995836, 0.06215885, 1.99436774, 0.77589616, 0.05003577, 2.45220238,
                1.14352085, 0.91552118, 0.43143425, 0.91820335, 0.33512254, 1.63188178, 0.45073285, 0.11618706, 0.98611168,
                0.6218516, 0.13076004, 0.49089057, 0.4967578, 0.63775862, 1.16009522, 0.29896259, 0.58792701, 1.153215,
                1.56923019, 0.54900448, 0.7352573, 0.87693621, 1.33628879, 0.44052683]
i = 0
j = 0
it = i

S_con = 0
trimen = 0

R[0] += 1
V[0] += arrive_time[i]

J5 = 0
JF5 = 0
Z5 = 0
Tl = 0

Ttime.append(arrive_time[i])

jN[it] = i + 1
jP[it] = arrive_time[i]
jQ[it] = 0
jQt[it] = 0
jS[it] = arrive_time[i]
jD[it] = service_time[0]
jF[it] = arrive_time[i] + service_time[0]
jk[it] = it + 1

Ttype.append(1)
condition.append(1)

# T1 += service_time[i]

Tremained.append(service_time[j])
Tnew.append(arrive_time[i + 1])
k[0] = service_time[j]
kj[0] = it+1
numj.append(it + 1)
numk.append(it + 1)
trimen = arrive_time[it + 1] - service_time[j]
S_con = 1
# V[S_con] += min(service_time[j],arrive_time[i + 1])
i += 1
it += 1
J5 += 1

while len(Ttime) != 100:
    if S_con == 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + Tnew[-1]
        jQ[i] = 0
        jQt[i] = 0
        jS[i] = Ttime[-1] + Tnew[-1]
        jk[i] = 1

        V[S_con] += Tnew[-1]
        Ttime.append(Ttime[-1] + Tnew[-1])
        Ttype.append(1)
        condition.append(S_con + 1)

        # T1 += service_time[i]

        Tremained.append(service_time[i])
        Tnew.append(arrive_time[i + 1])
        numj.append(i + 1)
        numk.append(1)
        trimen = arrive_time[i + 1] - service_time[j]
        S_con += 1
        R[S_con] += 1
        k[0] = service_time[j]
        kj[0] = i+1
        # V[S_con] += min(service_time[j], arrive_time[i + 1])
        it += 1
        i += 1
        J5 += 1
        Z5 += S_con
    elif trimen > 0:

        jS[kj[np.argmin(k)] - 1] = Ttime[-1] + k[np.argmin(k)] - service_time[kj[np.argmin(k)] - 1]
        jD[kj[np.argmin(k)] - 1] = service_time[kj[np.argmin(k)] - 1]
        jF[kj[np.argmin(k)] - 1] = Ttime[-1] + Tremained[-1]
        jk[kj[np.argmin(k)] - 1] = np.argmin(k) + 1

        V[S_con] += k[np.argmin(k)]
        numj.append(kj[np.argmin(k)])
        Ttime.append(Ttime[-1] + k[np.argmin(k)])
        k = np.subtract(k, k[np.argmin(k)])
        Ttype.append(2)
        condition.append(S_con - 1)
        j += 1
        if S_con <= n:
            kj[np.argmin(k)] = -1
            k[np.argmin(k)] = float('inf')
        else:
            k[np.argmin(k)] = float('inf')
            iop = np.argmax(k)
            jS[it] = Ttime[-1]
            jD[it] = service_time[it]
            jF[it] = Ttime[-1] + service_time[it]
            jk[it] = np.argmax(k) + 1
            k41[iop] += 1
            k42[iop] += service_time[kj[iop] - 1]
            k[iop] = service_time[it]
            kj[iop] = it + 1
            it+=1
            if it>i:
                i+=1
        if S_con - 1 == 0:
            Tremained.append(-1)
            Tnew.append(trimen)
            Tl += trimen
        else:
            Tremained.append(min(k))
            Tnew.append(trimen)

        trimen = trimen - min(k)
        numk.append(np.argmin(k) + 1)
        S_con -= 1
        R[S_con] += 1
        JF5 += 1
        Z5 += S_con
    elif trimen < 0:
        tp = trimen + min(k)
        k = np.subtract(k, tp)
        V[S_con] += tp
        jN[i] = i + 1
        jP[i] = Ttime[-1] + tp
        jQ[i] = S_con
        if S_con < n:
            jQt[i] = 0
        else:
            for d in range(S_con):
                jQt[i] += service_time[d + j + 1]
            jQt[i] -= trimen

        if max(k) == float('inf'):
            numk.append(np.argmax(k) + 1)
            jk[i] = np.argmax(k) + 1
            jS[it] = Ttime[-1] + tp + jQt[i]
            jD[it] = service_time[i]
            jF[it] = Ttime[-1] + tp + service_time[i]
            k41[np.argmax(k)] += 1
            k42[np.argmax(k)] += service_time[kj[np.argmax(k)] - 1]
            kj[np.argmax(k)] = i+1
            k[np.argmax(k)] = service_time[i]
            it = i+1
        else:
            numk.append(-1)

        Ttime.append(Ttime[-1] + tp)
        Ttype.append(1)
        condition.append(S_con + 1)
        Tremained.append(min(k))
        Tnew.append(arrive_time[i + 1])
        i += 1
        S_con += 1
        R[S_con] += 1
        numj.append(i)
        trimen = arrive_time[i] - min(k)
        J5 += 1
        Z5 += S_con

Rot = R / 100
print(sum(V))
print(Ttime[-1])
Vot = V / Ttime[-1]
Z5 = Z5 / 100
wer = 0
for d3 in range(len(jF)):
    if jF[d3] > 0:
        wer += (jF[d3] - jS[d3])
Tq5 = sum(jQt) / JF5
Tm5 = wer / JF5

for d1 in range(len(jQ)):
    if jQ[d1] < 5:
        jQ[d1] = 0
    else:
        jQ[d1] -= 4

for d2 in range(n):
    if k[d2] != float('inf'):
        k42[d2] -= k[d2]
    k43[d2] = Ttime[-1] - k42[d2]
    k44[d2] = k43[d2] / Ttime[-1]

f.write(str(np.around(L, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttime, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttype, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(condition, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tremained, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tnew, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(numj, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(numk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jN, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jP, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQ, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQt, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jD, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jF, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(R, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(V, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Rot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Vot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(R), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(V), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Rot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Vot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(J5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(JF5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Z5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tq5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Tm5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(k41, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k42, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k43, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k44, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')

k = [float('inf'), float('inf'), float('inf'), float('inf'), float('inf')]
kj = [-1, -1, -1, -1, -1]
m = 14
baned = []
Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []
numk = []

jN = np.zeros(100)
jP = np.zeros(100)
jQ = np.zeros(100)
jQt = np.zeros(100)
jS = np.zeros(100)
jD = np.zeros(100)
jF = np.zeros(100)
jk = np.zeros(100)

R = np.zeros(100)
V = np.zeros(100)

k41 = np.zeros(n)
k42 = np.zeros(n)
k43 = np.zeros(n)
k44 = np.zeros(n)

T0 = Tt
T1 = 0

arrive_time = expon.rvs(scale=1 / lambd, size=105)
print(arrive_time)
service_time = expon.rvs(scale=1 / u, size=105)
print(service_time)
i = 0
j = 0
it = i

S_con = 0
trimen = 0

R[0] += 1
V[0] += arrive_time[i]

J5 = 0
JF5 = 0
Z5 = 0
Tl = 0

Ttime.append(arrive_time[i])

jN[it] = i + 1
jP[it] = arrive_time[i]
jQ[it] = 0
jQt[it] = 0
jS[it] = arrive_time[i]
jD[it] = service_time[0]
jF[it] = arrive_time[i] + service_time[0]
jk[it] = it + 1

Ttype.append(1)
condition.append(1)

# T1 += service_time[i]

Tremained.append(service_time[j])
Tnew.append(arrive_time[i + 1])
k[0] = service_time[j]
kj[0] = it+1
numj.append(it + 1)
numk.append(it + 1)
trimen = arrive_time[it + 1] - service_time[j]
S_con = 1
# V[S_con] += min(service_time[j],arrive_time[i + 1])
i += 1
it += 1
J5 += 1

while len(Ttime) != 100:
    if S_con == 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + Tnew[-1]
        jQ[i] = 0
        jQt[i] = 0
        jS[i] = Ttime[-1] + Tnew[-1]
        jk[i] = 1

        V[S_con] += Tnew[-1]
        Ttime.append(Ttime[-1] + Tnew[-1])
        Ttype.append(1)
        condition.append(S_con + 1)

        # T1 += service_time[i]

        Tremained.append(service_time[j])
        Tnew.append(arrive_time[i + 1])
        numj.append(i + 1)
        numk.append(1)
        trimen = arrive_time[i + 1] - service_time[j]
        S_con += 1
        R[S_con] += 1
        k[0] = service_time[j]
        kj[0] = i+1
        it += 1
        while it in baned:
            it += 1
        i += 1
        J5 += 1
        Z5 += S_con
    elif trimen > 0:

        jS[kj[np.argmin(k)] - 1] = Ttime[-1] + k[np.argmin(k)] - service_time[kj[np.argmin(k)] - 1]
        jD[kj[np.argmin(k)] - 1] = service_time[kj[np.argmin(k)] - 1]
        jF[kj[np.argmin(k)] - 1] = Ttime[-1] + Tremained[-1]
        jk[kj[np.argmin(k)] - 1] = np.argmin(k) + 1

        V[S_con] += k[np.argmin(k)]
        numj.append(kj[np.argmin(k)])
        Ttime.append(Ttime[-1] + k[np.argmin(k)])
        k = np.subtract(k, k[np.argmin(k)])
        Ttype.append(3)
        condition.append(S_con - 1)
        j += 1
        if S_con <= n:
            kj[np.argmin(k)] = -1
            k[np.argmin(k)] = float('inf')
        else:
            k[np.argmin(k)] = float('inf')
            iop = np.argmax(k)
            k41[iop] += 1
            k42[iop] += service_time[kj[iop] - 1]
            k[iop] = service_time[it]
            kj[iop] = it + 1
            it += 1
            while it in baned:
                it += 1
            if it > i:
                i += 1
        if S_con - 1 == 0:
            Tremained.append(-1)
            Tnew.append(trimen)
            Tl += trimen
        else:
            Tremained.append(min(k))
            Tnew.append(trimen)

        trimen = trimen - min(k)
        numk.append(np.argmin(k) + 1)
        S_con -= 1
        R[S_con] += 1
        JF5 += 1
        Z5 += S_con
    elif trimen < 0:
        if S_con < n + m:
            tp = trimen + min(k)
            k = np.subtract(k, tp)
            V[S_con] += tp
            jN[i] = i + 1
            jP[i] = Ttime[-1] + tp
            jQ[i] = S_con
            if S_con < n:
                jQt[i] = 0
            else:
                for d in range(S_con):
                    jQt[i] += service_time[d + j + 1]
                jQt[i] -= trimen

            if max(k) == float('inf'):
                numk.append(np.argmax(k) + 1)
                jk[i] = np.argmax(k) + 1
                k41[np.argmax(k)] += 1
                k42[np.argmax(k)] += service_time[kj[np.argmax(k)] - 1]
                kj[np.argmax(k)] = i + 1
                k[np.argmax(k)] = service_time[i]
                it = i+1
            else:
                numk.append(-1)

            Ttime.append(Ttime[-1] + tp)
            Ttype.append(1)
            condition.append(S_con + 1)
            Tremained.append(min(k))
            Tnew.append(arrive_time[i + 1])
            i += 1
            S_con += 1
            R[S_con] += 1
            numj.append(i)
            trimen = arrive_time[i] - min(k)
            J5 += 1
            Z5 += S_con
        else:
            tp = trimen + min(k)
            k = np.subtract(k, tp)
            V[S_con] += tp
            jN[i] = i + 1
            jP[i] = Ttime[-1] + tp
            jQ[i] = -1
            jQt[i] = 0
            jS[it] = -1
            jD[it] = 0
            jF[it] = jP[i]
            jk[it] = -2

            Ttime.append(Ttime[-1] + tp)
            Ttype.append(2)
            numk.append(-1)
            condition.append(S_con)
            Tremained.append(min(k))
            Tnew.append(arrive_time[i + 1])
            baned.append(i)
            i += 1
            R[S_con] += 1
            # V[S_con] += min(-trimen, arrive_time[i])
            numj.append(i)
            trimen = arrive_time[i] - min(k)
            J5 += 1
            Z5 += S_con

Rot = R / 100
print(sum(V))
print(Ttime[-1])
Vot = V / Ttime[-1]
Z5 = Z5 / 100
wer = 0
for d3 in range(len(jF)):
    if jF[d3] > 0:
        wer += (jF[d3] - jS[d3])
Tq5 = sum(jQt) / JF5
Tm5 = wer / JF5

for d1 in range(len(jQ)):
    if jQ[d1] < 5:
        jQ[d1] = 0
    else:
        jQ[d1] -= 4

for d2 in range(n):
    if k[d2] != float('inf'):
        k42[d2] -= k[d2]
    k43[d2] = Ttime[-1] - k42[d2]
    k44[d2] = k43[d2] / Ttime[-1]

f.write(str(np.around(L, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttime, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Ttype, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(condition, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tremained, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tnew, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(numj, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(numk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jN, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jP, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQ, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jQt, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jD, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jF, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(jk, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(R, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(V, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Rot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Vot, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(R), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(V), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Rot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(sum(Vot), 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(J5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(JF5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Z5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tq5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(Tm5, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(k41, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k42, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k43, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(k44, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')

ro = lambd / u
v = ro / m
v1 = ro/n
print('ro')
print(ro)
Ro1 = 1/(1 + ro + ro * ro / 2 + ro * ro * ro / 6 + ro * ro * ro * ro / 24 + ((ro * ro * ro * ro * ro / 120) * (1/(1-v1))))
print('Ro')
print(Ro1)
Ro2 = 1/(1 + ro + ro * ro / 2 + ro * ro * ro / 6 + ro * ro * ro * ro / 24 + (ro * ro * ro * ro * ro / 120) * (
            1 + v + pow(v, 2) + pow(v, 3) + pow(v, 4) + pow(v, 5)
            + pow(v, 6) + pow(v, 7) + pow(v, 8) + pow(v, 9)
            + pow(v, 10) + pow(v, 11) + pow(v, 12) + pow(v, 13)
            + pow(v, 14)))
K1 = ro
K2 = ro * (1 - ((pow(ro, 5) / 120) * pow(v, m-n) * Ro2))
print(Ro2)
Q1 = v * pow(ro, n) / 120 * Ro1 / ((1 - v) * (1 - v))
Q2 = v * pow(ro, n) / 120 * Ro2 * (1 - (m + 1) * pow(v, m) + m * pow(v, m + 1)) / ((1 - v) * (1 - v))
Z1 = K1 + Q1
Z2 = K2 + Q2
tq1 = Q1 / lambd
tq2 = Q2 / lambd
tin1 = Z1 / lambd
tin2 = Z2 / lambd



f.write(str(np.around(K1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Q1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Z1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(tq1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(tin1, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(K2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Q2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Z2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(tq2, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(tin2, 5)))
f.write('\n')
f.write('\n')
t = Ro2 * pow(ro, n) * pow(v, m) / 120
print(t)
f.write(str(np.around(t, 16)))
f.write('\n')
f.write(str(np.around(t, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
r5 = Ro1*pow(ro,5)/math.factorial(5)
l = 0
for i in range(6):
    print(i)
    print(Ro1*pow(ro,i)/math.factorial(i))
    l += (Ro1*pow(ro,i)/math.factorial(i))

print('r5')
print(r5)
for i in range(6,13):
    print(i)
    print(r5*pow(ro/n,i-n))
    l += (r5*pow(ro/n,i-n))
print(l)