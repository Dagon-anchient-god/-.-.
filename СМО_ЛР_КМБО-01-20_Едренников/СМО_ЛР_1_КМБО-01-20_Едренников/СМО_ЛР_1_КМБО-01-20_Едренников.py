import numpy as np
from scipy.stats import expon
from decimal import Decimal

f = open('answer.txt', 'r+')

u = 1.254
lambd = 1.052
Ts = 0.794
Tt = 0.784  # Время между заявками

L = list(range(1, 101))
Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []

jN = np.zeros(100)
jP = np.zeros(100)
jS = np.zeros(100)
jF = np.zeros(100)

T0 = Tt
T1 = 0

J5 = 0
JF5 = 0
JL5 = 0
Dl = 0
Tl = 0

service_time = expon.rvs(scale=1 / u, size=100)
service_time = [0.0596276, 0.72127201, 2.01710223, 0.32372014, 0.44929662, 1.38898635, 0.28999866, 2.82485371,
                0.00489021, 0.80920563, 0.53380782, 0.33801901, 0.82054915, 0.25590868, 1.18827375, 0.26642167, 1.7424584,
                0.07252952, 0.212955, 1.53165004, 0.00627624, 1.17975602, 1.00696473, 1.02885748, 1.1812449, 0.99447012,
                1.0626488, 0.85383117, 0.71700256, 0.3797618, 0.02984504, 0.36360791, 0.16964442, 0.49890133, 0.0947652,
                1.22494775, 4.03252267, 0.14152315, 0.50042083, 0.39265258, 1.88307502, 0.2069792, 0.85826512, 0.80996065,
                1.12609482, 1.31039982, 1.51896547, 0.85907068, 1.14635457, 1.60915022, 0.17615997, 0.15130618, 0.15950662,
                0.0796058, 0.35177219, 0.18433293, 0.6734803, 1.12783108, 0.04402474, 1.18628942, 0.02226666, 0.1077807, 1.06877713,
                0.0984029, 0.27811944, 0.51479237, 0.38918731, 0.13941547, 0.19589416, 0.35329383, 0.0173077, 1.53808859,
                0.023749, 0.10802429, 0.09758635, 0.38637948, 1.23661508, 1.059142, 0.14084473, 0.01654323, 1.26470771, 0.28802272,
                0.03494353, 0.33704994, 0.16819181, 0.5593081, 0.29836402, 0.46637185, 0.97799658, 0.44816637, 0.16128861, 1.01136497,
                2.1438003, 0.16640893, 0.15373774, 0.02907984, 1.05664218, 0.32663872, 1.32075165, 1.1589173]
print(service_time)
i = 0
j = 0

jN[j] = j + 1
jP[j] = Tt
jS[j] = service_time[i]
Ttime.append(Tt)



Ttype.append(1)
condition.append(1)

T1 += service_time[i]

Tremained.append(service_time[i])
Tnew.append(Tt)
j += 1
jr = j
numj.append(j)
trimen = service_time[i] - Tt
S_con = 1
J5+=1
while len(Ttime) != 100:
    if S_con == 0:
        T0 += -trimen

        jN[j] = j + 1
        jP[j] = Ttime[-1] - trimen
        jS[j] = service_time[i]

        Ttime.append(Ttime[-1] - trimen)
        Ttype.append(1)
        condition.append(1)

        T1 += service_time[i]

        Tremained.append(service_time[i])
        Tnew.append(Tt)
        j += 1
        jr = j
        numj.append(j)
        trimen = service_time[i] - Tt
        S_con = 1
        J5 += 1
    elif trimen < 0 and Ttype[-1] == 2:

        jF[jr-1] = Ttime[-1] + Tremained[-1]

        Ttime.append(Ttime[-1] + Tremained[-1])
        Ttype.append(3)
        condition.append(0)
        Tremained.append(-1)

        T0 += -trimen

        Tnew.append(-trimen)
        Tl += -trimen
        numj.append(jr)
        i += 1
        S_con = 0
        JF5 += 1
    elif S_con == 1 and trimen > 0:

        jN[j] = j + 1
        jP[j] = Ttime[-1] + Tt
        jS[j] = 0
        jF[j] = Ttime[-1] + Tt

        Ttime.append(Ttime[-1] + Tt)
        Ttype.append(2)
        condition.append(1)
        Tremained.append(trimen)
        Tnew.append(Tt)
        j += 1
        numj.append(j)
        trimen -= Tt
        J5 += 1
        JL5 += 1
    elif trimen < 0:

        jF[jr-1] = Ttime[-1] + Tremained[-1]

        Ttime.append(Ttime[-1] + Tremained[-1])
        Ttype.append(3)
        condition.append(0)
        Tremained.append(-1)

        T0 += -trimen

        Tnew.append(-trimen)
        numj.append(jr)
        i += 1
        S_con = 0
        Tl += -trimen
        JF5 += 1

R0 = condition.count(0)
R1 = condition.count(1)
v0 = R0 / 100
v1 = R1 / 100
T0 = Ttime[-1] - T1
delta0 = T0 / Ttime[-1]
delta1 = T1 / Ttime[-1]
Dl = JL5/J5
Tl = Tl/Ttime[-1]

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
f.write(str(np.around(jN, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jP, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jS, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(jF, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(R0, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(R1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(v0, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(v1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(T0, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(T1, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(delta0, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(delta1, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(np.around(J5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(JF5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(JL5, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Dl, 5)))
f.write('\n')
f.write('\n')
f.write(str(np.around(Tl, 5)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')
f.write('\n')

Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []

jN = np.zeros(100)
jP = np.zeros(100)
jQ = np.zeros(100)
jQt = np.zeros(100)
jS = np.zeros(100)
jD = np.zeros(100)
jF = np.zeros(100)

R = np.zeros(100)
V = np.zeros(100)

T0 = Tt
T1 = 0

arrive_time = expon.rvs(scale=1 / lambd, size=100)
arrive_time = [0.40117091, 1.0135882, 3.16533116, 0.03230948, 1.26603334, 0.80838107, 0.01125001, 0.69512731, 0.68403592, 0.90059743, 1.65994806, 1.24489832, 0.17155369, 0.28960933, 1.70368105, 1.04902514, 1.56596271, 0.18996548, 1.41560046, 1.01285334, 0.07614671, 2.87264663, 0.3635338, 0.9268059, 0.30303236, 0.0857338, 0.06229516, 0.7802278, 1.52775425, 0.12170897, 0.60564453, 0.55114783, 1.64119539, 0.6147427, 0.08273187, 0.37517914, 0.00815103, 3.16347353, 1.23067033, 0.24756659, 0.17842878, 0.43756393, 0.72285051, 1.18334614, 1.02804347, 1.86880286, 0.68978814, 0.88748695, 0.10272516, 0.74443909, 0.7709879, 1.05304843, 0.74895205, 0.88201039, 1.55375056, 0.02729069, 0.89429496, 0.65826222, 0.25342226, 0.25016314, 1.14511598, 1.14859708, 0.94699604, 0.63564956, 0.44301601, 0.18424226, 0.48948496, 0.21349709, 1.06094758, 0.14500643, 0.99892129, 3.42469012, 0.85917016, 2.33556292, 0.72154634, 0.18810993, 0.36188151, 0.13434438, 1.99659399, 1.41253129, 1.26997202, 0.28451984, 0.93023139, 0.62328968, 0.89739204, 1.07336679, 0.03165158, 1.46436741, 0.48715263, 1.61223369, 0.2703513, 0.34097123, 0.38275818, 1.89983983, 2.47692569, 0.60724232, 0.83825833, 1.27701367, 0.37523242, 0.11070942]
print(arrive_time)
i = 0
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

Ttype.append(1)
condition.append(1)

# T1 += service_time[i]

Tremained.append(Ts)
Tnew.append(arrive_time[i + 1])
numj.append(it + 1)
trimen = arrive_time[it + 1] - Ts
S_con = 1
R[S_con] += 1
#V[S_con] += min(Ts,arrive_time[i + 1])
i+=1
it+=1
J5 += 1
while len(Ttime) != 100:
    if S_con == 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + Tnew[-1]
        jQ[i] = 0
        jQt[i] = 0
        jS[i] = Ttime[-1] + Tnew[-1]

        V[S_con] += Tnew[-1]
        Ttime.append(Ttime[-1] + Tnew[-1])
        Ttype.append(1)
        condition.append(S_con + 1)

        # T1 += service_time[i]

        Tremained.append(Ts)
        Tnew.append(arrive_time[it + 1])
        numj.append(it + 1)
        trimen = arrive_time[it + 1] - Ts
        S_con += 1
        R[S_con] += 1
        #V[S_con] += min(Ts, arrive_time[it + 1])
        it += 1
        i += 1
        J5 += 1
        Z5 += S_con
    elif trimen > 0:

        jS[it] = Ttime[-1] + Tnew[-1] - Ts
        jD[it] = Ts
        jF[it] = Ttime[-1] + Tremained[-1]

        #V[S_con] += Tremained[-1]
        Ttime.append(Ttime[-1] + Tremained[-1])
        Ttype.append(2)
        V[S_con] += Tremained[-1]
        condition.append(S_con - 1)
        if S_con - 1 == 0:
            Tremained.append(-1)
            Tnew.append(trimen)
            #V[S_con-1] += trimen
            Tl += trimen
        else:
            #V[S_con-1] += Tremained[-1]
            Tremained.append(Ts)
            Tnew.append(0)
            #V[S_con-1] += min(Ts, trimen)

        # T0 += -trimen

        trimen = arrive_time[i + 1] - Ts
        numj.append(it)
        if S_con - 1 != 0:
            it += 1
        S_con -= 1
        R[S_con] += 1
        JF5 += 1
        Z5 += S_con
    elif trimen < 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + arrive_time[i]
        jQ[i] = S_con
        jQt[i] = Ts * (S_con) + trimen

        V[S_con] += arrive_time[i]
        Ttime.append(Ttime[-1] + arrive_time[i])
        Ttype.append(1)
        condition.append(S_con + 1)
        Tremained.append(-trimen)
        Tnew.append(arrive_time[i + 1])
        i += 1
        S_con += 1
        R[S_con] += 1
        #V[S_con] += min(-trimen, arrive_time[i])
        numj.append(i)
        trimen += arrive_time[i]
        J5 += 1
        Z5 += S_con





Rot = R/100
Vot = V/Ttime[-1]
print(sum(V))
print(Ttime[-1])
Z5 = Z5/100
Tq5 = sum(jQt)/JF5
Tm5 = sum(jF-jP)/JF5
Tl = Tl/Ttime[-1]

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
f.write(str(np.around(Tl, 5)))
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
print(sum(V))



Ttime = []
Ttype = []
condition = []
Tremained = []
Tnew = []
numj = []

jN = np.zeros(100)
jP = np.zeros(100)
jQ = np.zeros(100)
jQt = np.zeros(100)
jS = np.zeros(100)
jD = np.zeros(100)
jF = np.zeros(100)

R = np.zeros(100)
V = np.zeros(100)

T0 = Tt
T1 = 0

arrive_time = expon.rvs(scale=1 / lambd, size=105)
arrive_time = [2.20718492, 0.26529921, 0.10325718, 0.12288507, 2.39239022, 0.22729359, 0.71908768, 0.19363402, 1.31226468, 0.59389162, 0.24181767, 0.02464913, 0.45277553, 0.60847666, 0.02775334, 1.13163481, 1.5844468, 3.03362209, 0.43993822, 0.23238649, 0.43077442, 0.4603215, 0.08828852, 0.16520653, 0.31596391, 1.37001638, 0.56499403, 0.43320493, 0.38779847, 0.41418553, 3.09062796, 1.00399573, 1.16567395, 1.6422182, 1.05017879, 0.17518421, 2.14681417, 0.23755632, 2.56018167, 0.15162922, 0.3140315, 0.19054134, 0.11435043, 0.53299786, 2.18066243, 0.80918078, 0.60017746, 1.70391298, 4.07399, 3.30689993, 1.02838463, 0.31783687, 1.03291466, 0.11887342, 0.14481511, 0.07579462, 0.63721745, 1.95042728, 0.15373085, 2.49893149, 1.51821568, 2.30314229, 0.31068136, 1.45814779, 0.80388062, 0.60584391, 0.89974687, 0.06774533, 0.5994367, 0.15957083, 0.82415271, 0.35558941, 0.54138821, 3.28102585, 0.23260546, 2.12848419, 0.99654792, 0.33818901, 0.11495295, 0.04273647, 0.03962951, 0.17755221, 0.61092236, 2.19789821, 1.02850586, 0.44404132, 0.80435907, 1.84632329, 2.56700775, 2.23089819, 0.2723861, 0.97549922, 0.34930668, 0.18054849, 0.69297077, 0.79819385, 4.04333174, 0.3659875, 0.91265219, 2.01531124, 0.05463354, 0.45549287, 1.52315717, 1.2002214, 0.7876374,]
print(arrive_time)
service_time = expon.rvs(scale=1 / u, size=105)
service_time = [2.93613062e+00, 1.56938416e+00, 1.52941779e+00, 8.39128275e-02, 5.91677518e-01, 3.55420398e-02, 1.44913065e+00, 5.53054285e-01, 3.01199278e+00, 2.35089535e+00, 3.84852334e-01, 5.04537144e-01, 9.46109388e-01, 1.26553492e+00, 4.93639929e-01, 3.90237738e-01, 3.37337979e+00, 3.11337732e+00, 8.85186315e-01, 5.87309448e-02, 4.12249674e-01, 2.21708319e-01, 1.12858638e+00, 9.02687314e-02, 1.44916871e-01, 5.65828381e-02, 7.11930058e-01, 9.04080382e-01, 4.02732241e+00, 4.37148373e-01, 1.60947711e-01, 2.61202992e-01, 3.03890393e-02, 5.64698266e-01, 1.00128821e+00, 2.72552891e-01, 1.28101678e+00, 1.77751875e-01, 7.50633058e-02, 1.15748907e-01, 1.57348018e+00, 3.44341645e-01, 8.37761906e-01, 2.35317669e-01, 1.32429834e-01, 3.51273978e-01, 9.83436506e-01, 2.87693584e-01, 1.69060553e+00, 2.12867780e+00, 1.74323242e+00, 9.39580009e-01, 2.01175799e-01, 5.71501694e-02, 4.47566636e-01, 2.21465944e+00, 1.25437703e+00, 2.82898191e-01, 1.55794780e+00, 5.25302804e-01, 1.49683128e+00, 8.78322611e-01, 6.94661885e-02, 1.30313830e-01, 3.02726805e-01, 3.19998528e+00, 4.81134657e-01, 7.21457649e-01, 4.93482834e-02, 4.36012464e-01, 6.39453990e-01, 8.99441309e-02, 2.45257555e+00, 1.70642016e+00, 6.14330704e-01, 9.13332468e-01, 2.26219077e-01, 1.93888673e-01, 1.30890614e+00, 1.38956501e-03, 5.43608367e-01, 6.39416687e-02, 9.02842252e-01, 2.11538645e+00, 3.80504429e-02, 2.69041841e-01, 1.20720293e+00, 5.36911339e-01, 1.17333053e+00, 1.34255053e-01, 8.01150684e-01, 8.24543871e-01, 8.39682266e-01, 1.32277362e+00, 1.95569786e-01, 1.08849715e-01, 2.01899368e-01, 1.73524631e+00, 3.05655782e-01, 2.76517344e-01, 5.15472337e-01, 9.14255974e-01, 1.54729761e+00, 2.49196005e-01, 1.50307918e+00,]
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

Ttype.append(1)
condition.append(1)

# T1 += service_time[i]

Tremained.append(service_time[j])
Tnew.append(arrive_time[i + 1])
numj.append(it + 1)
trimen = arrive_time[it + 1] - service_time[j]
S_con = 1
R[S_con] += 1
#V[S_con] += min(service_time[j],arrive_time[i + 1])
i+=1
it+=1
J5 += 1
while len(Ttime) != 100:
    if S_con == 0:

        jN[i] = i + 1
        jP[i] = Ttime[-1] + Tnew[-1]
        jQ[i] = 0
        jQt[i] = 0
        jS[i] = Ttime[-1] + Tnew[-1]

        V[S_con] += Tnew[-1]
        Ttime.append(Ttime[-1] + Tnew[-1])
        Ttype.append(1)
        condition.append(S_con + 1)

        # T1 += service_time[i]

        Tremained.append(service_time[j])
        Tnew.append(arrive_time[it + 1])
        numj.append(it + 1)
        trimen = arrive_time[it + 1] - service_time[j]
        S_con += 1
        R[S_con] += 1
        #V[S_con] += min(service_time[j], arrive_time[i + 1])
        it += 1
        i += 1
        J5 += 1
        Z5 += S_con
    elif trimen > 0:

        jS[it] = Ttime[-1] + Tnew[-1] - service_time[j]
        jD[it] = service_time[j]
        jF[it] = Ttime[-1] + Tremained[-1]

        Ttime.append(Ttime[-1] + Tremained[-1])
        Ttype.append(2)
        V[S_con] += Tremained[-1]
        condition.append(S_con - 1)
        j += 1
        if S_con - 1 == 0:
            Tremained.append(-1)
            Tnew.append(trimen)
            #V[S_con-1] += trimen
            Tl += trimen
        else:
            #V[S_con - 1] += Tremained[-1]
            Tremained.append(service_time[j])
            Tnew.append(0)
            #V[S_con-1] += min(service_time[j], trimen)

        # T0 += -trimen

        trimen = arrive_time[i + 1] - service_time[j]
        numj.append(it)
        if S_con - 1 != 0:
            it += 1
        S_con -= 1
        R[S_con] += 1
        JF5 += 1
        Z5 += S_con
    elif trimen < 0:

        V[S_con] += arrive_time[i]
        jN[i] = i + 1
        jP[i] = Ttime[-1] + arrive_time[i]
        jQ[i] = S_con
        for k in range(S_con):
            jQt[i] += service_time[j+k+1]
        jQt[i] += trimen

        Ttime.append(Ttime[-1] + arrive_time[i])
        Ttype.append(1)
        condition.append(S_con + 1)
        Tremained.append(-trimen)
        Tnew.append(arrive_time[i + 1])
        i += 1
        S_con += 1
        R[S_con] += 1
        #V[S_con] += min(-trimen, arrive_time[i])
        numj.append(i)
        trimen += arrive_time[i]
        J5 += 1
        Z5 += S_con



Rot = R/100
print(sum(V))
print(Ttime[-1])
Vot = V/Ttime[-1]
Z5 = Z5/100
Tq5 = sum(jQt)/JF5
Tm5 = sum(jF-jP)/JF5
Tl = Tl/Ttime[-1]


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
f.write(str(np.around(Tl, 5)))
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
print(sum(V))