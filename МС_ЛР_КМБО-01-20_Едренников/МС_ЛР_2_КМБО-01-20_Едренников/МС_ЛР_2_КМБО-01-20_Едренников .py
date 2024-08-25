import numpy as np
import math
import scipy.stats
from scipy.stats import hypergeom
from scipy.stats import binom
from scipy.stats import norm
from scipy.stats import expon
from scipy.stats import uniform
from scipy.stats import chi2
import statistics
import matplotlib.pyplot as plt
import pandas as pd

f = open('answer.txt', 'r+')
a = 0.60
b = 3.06
data_unif = [0.87237, 0.61947, 2.28114, 1.05954, 3.01849, 0.78614, 1.97266, 0.81233, 0.66475, 1.17367, 0.81927, 0.92895,
             2.00952, 0.92135, 1.97646, 1.67131, 1.32465, 2.54559, 0.87393, 1.66786,
             2.86762, 0.88431, 0.91213, 3.04601, 3.00578, 1.85865, 2.45107, 1.17334, 0.73715, 2.23896, 1.34590, 1.21325,
             2.25551, 1.23435, 3.01429, 1.13780, 1.76009, 1.09206, 1.35849, 2.43128,
             2.71963, 2.29392, 2.72527, 2.91515, 1.44885, 2.62008, 1.04418, 1.30352, 2.46744, 0.89681, 1.15615, 2.02147,
             1.07346, 1.74173, 2.86102, 2.32731, 2.71160, 2.51601, 2.23878, 1.46360,
             2.54880, 1.85966, 1.14817, 2.64072, 2.99352, 1.50133, 1.01462, 2.11928, 1.81440, 1.41270, 1.14874, 1.98375,
             1.49910, 2.69859, 2.99160, 2.17237, 1.07898, 2.73331, 1.92270, 2.66263,
             1.81616, 2.58813, 0.60401, 1.79972, 0.76114, 2.55662, 1.52095, 1.32455, 0.80442, 1.52681, 1.22406, 2.45410,
             0.86988, 2.50577, 0.89397, 1.50328, 1.54239, 2.73251, 2.14441, 1.65987,
             1.53932, 2.25165, 1.39593, 2.01633, 1.40785, 0.77011, 1.94263, 2.21776, 1.41710, 2.98715, 2.00098, 1.74565,
             2.78562, 0.70137, 1.72870, 1.19845, 0.97294, 1.43609, 0.79789, 1.69975,
             2.97624, 1.78702, 0.96663, 1.89037, 2.91193, 1.20819, 2.94638, 1.97961, 1.37576, 2.36926, 2.84769, 1.63859,
             0.96078, 1.75007, 2.34675, 1.55341, 2.82976, 2.09428, 1.40535, 2.90010,
             0.97956, 2.26028, 2.43269, 1.48038, 2.06027, 0.80024, 2.20077, 2.85829, 1.04978, 2.40315, 1.88943, 1.62652,
             0.89097, 1.37580, 1.91245, 0.76138, 1.25136, 2.08915, 1.36937, 1.93149,
             1.13584, 1.90767, 2.50076, 1.02003, 2.18454, 2.18617, 2.37443, 2.43942, 0.99391, 2.66653, 1.07565, 1.76319,
             2.44607, 1.81466, 1.90012, 1.26679, 1.62841, 1.41934, 1.19940, 0.87180,
             0.95918, 2.64926, 2.72230, 0.69900, 1.52368, 0.94387, 0.87747, 2.40600, 0.83509, 2.09499, 1.77365, 2.11390,
             2.47342, 1.47479, 3.05531, 1.14694, 1.59200, 1.85078, 1.88687, 2.05778]
data_unif.sort()
data_exp = [0.97494, 0.54345, 0.55448, 0.89675, 0.16532, 0.14625, 0.27871, 0.08319, 0.57546, 0.75401,
            0.33935, 0.07207, 0.86823, 2.27584, 1.43015, 0.82341, 2.02619, 0.29441, 0.91136, 0.65958,
            1.66464, 0.87931, 0.28180, 0.99399, 0.24299, 0.21930, 1.56621, 0.07542, 0.44137, 0.90047,
            0.61899, 1.97023, 0.18752, 0.28146, 0.61500, 0.30781, 0.21578, 0.23441, 0.23171, 0.49810,
            0.98031, 0.35218, 0.58011, 1.08712, 0.14496, 0.15921, 0.39636, 0.28065, 1.83692, 0.06197,
            0.33707, 0.88832, 0.39735, 0.50391, 0.44614, 0.65537, 0.98366, 0.16606, 1.86181, 0.94307,
            0.49935, 2.21054, 3.47755, 0.64380, 0.01571, 0.48544, 0.22525, 0.02169, 0.09666, 0.00583,
            0.26713, 1.20799, 0.18829, 0.14729, 0.06153, 0.60877, 2.18966, 0.31837, 0.26445, 0.48635,
            0.01271, 1.24150, 0.61427, 0.96416, 0.32586, 0.02224, 1.14577, 0.40663, 0.28343, 0.46425,
            0.67296, 0.59716, 0.36653, 0.88578, 0.16422, 0.03676, 0.44591, 0.47940, 0.16182, 0.57747,
            0.24873, 0.05087, 0.45316, 0.19181, 0.76124, 1.46750, 0.84597, 0.76690, 0.81015, 0.70416,
            0.87103, 1.43720, 0.16612, 0.55884, 2.15408, 0.08889, 0.32100, 0.54827, 0.19218, 0.31592,
            0.94434, 0.23889, 0.33519, 0.40694, 0.05246, 1.54816, 0.79350, 1.76549, 0.42547, 0.39379,
            0.15872, 0.22889, 1.46276, 1.04608, 0.76953, 0.11451, 0.59191, 0.31786, 1.43878, 1.05002,
            2.05666, 0.66660, 0.25338, 0.33463, 0.67874, 0.94867, 0.08803, 1.21272, 0.11179, 2.37156,
            0.51527, 0.12643, 0.02899, 1.92154, 1.35496, 0.55238, 0.20668, 0.01712, 1.06065, 0.88205,
            0.10873, 1.26062, 0.49399, 0.80641, 0.25320, 1.24723, 2.37736, 0.74403, 0.63150, 0.73679,
            0.47896, 0.64002, 0.32848, 0.65269, 0.70420, 0.47792, 1.36441, 0.15506, 0.07258, 0.22047,
            0.24231, 2.08651, 0.90103, 1.29706, 0.21108, 0.01437, 0.48383, 0.19103, 0.02780, 1.56642,
            0.03986, 0.75419, 0.68418, 1.27559, 0.11814, 0.02245, 0.17843, 1.68521, 0.72983, 0.76955]
data_exp.sort()
data_norm = [2.86728, 0.64330, 0.56871, -0.51954, 1.88552, 1.87000, 2.54401, 0.77150, -0.25163, 0.45270,
             1.13188, -0.13516, 2.79170, 1.54254, 0.86840, -1.67046, 0.24592, 1.38639, 0.26333, 0.54415,
             1.65331, 0.29311, -0.51887, 2.48220, 1.01074, 1.96870, 3.25616, 0.99819, 0.83646, 3.19009,
             -0.24485, -0.57104, 1.53098, 0.24618, 2.25014, 1.52737, 1.48454, 1.45196, 0.49125, 1.66006,
             0.77905, 0.11235, 1.39604, 3.44068, 0.21113, 2.91672, -0.41808, 2.38678, 2.40915, 2.38930,
             2.30288, 1.26531, -0.25179, 1.40918, 3.38322, 0.65500, 1.22189, 2.88522, 0.91645, -0.53520,
             2.63247, 0.46257, 2.48092, 2.70751, 0.26832, 1.89455, -0.32632, 1.05800, 2.39197, 1.45039,
             -0.68905, 2.35009, 2.74350, 2.04859, 2.18862, 1.85233, 0.99772, 0.78952, 0.87469, 0.99794,
             0.31365, 1.24173, -0.57746, 1.93329, 1.41590, 1.22422, 1.73337, 3.04611, 2.53595, -0.12383,
             3.83315, 1.51367, 4.21158, 4.26113, 3.74643, 0.35489, 2.17630, 2.74079, -0.50449, 0.83012,
             0.78376, 2.39267, 2.04421, 0.98281, 1.85254, 0.01681, 0.98511, 1.57608, 1.26465, 2.91597,
             0.53350, 2.97226, 1.64747, 2.77579, 1.89312, 1.52409, 2.13136, 0.70852, 0.86161, 1.64958,
             2.64152, 1.18987, 3.30493, 0.26305, 1.25149, 0.49265, 0.77461, 2.75295, 0.20638, 3.23153,
             -0.41601, 0.49294, 3.57286, 1.78484, 2.95677, 2.34132, 2.61648, 0.72660, 2.11193, 0.50248,
             2.47467, 1.34834, -0.24358, 0.35775, 2.69114, 0.27944, 1.78820, 1.54063, 2.77200, 1.66152,
             0.87717, 2.37526, 1.65645, 2.36921, 1.07254, 1.07622, 0.69913, -0.65297, 1.22832, 3.34510,
             -0.04131, 0.53450, 0.36705, 0.33717, 0.75314, 3.35553, -0.36351, 0.69532, 2.88261, -0.17766,
             1.26223, 1.61054, 1.76563, 2.04095, -0.41430, 0.70619, 0.41194, 2.50029, 1.30179, 0.86940,
             0.79233, -1.60799, 4.39980, 0.70509, 1.65828, -1.00514, 0.69058, 1.17355, 1.45649, 1.63657,
             1.58629, 1.50408, 2.46199, 0.30441, 1.59389, 3.32481, 0.78103, 3.11943, 2.11255, 3.04706]
data_norm.sort()
m = 1 + math.floor(math.log2(200))
a1 = [None] * (m + 1)
x1 = []

a2 = [None] * (m + 1)
x2 = []

a3 = [None] * (m + 1)
x3 = []

a1[0] = 0
a1[m] = max(data_exp)

a2[0] = min(data_norm)
a2[m] = max(data_norm)

a3[0] = a
a3[m] = b

for i in range(1, m):
    a1[i] = a1[i - 1] + (a1[m] - a1[0]) / m
    a2[i] = a2[i - 1] + (a2[m] - a2[0]) / m
    a3[i] = a3[i - 1] + (a3[m] - a3[0]) / m
a1.sort()
a2.sort()
a3.sort()
for i in range(1, m + 1):
    x1.append((a1[i - 1] + a1[i]) / 2)
    x2.append((a2[i - 1] + a2[i]) / 2)
    x3.append((a3[i - 1] + a3[i]) / 2)

pan = pd.Series(data_exp)
frequency1 = pan.groupby(pd.cut(pan, bins=a1, right=True)).count()
frequency1 = frequency1.tolist()
relative_frequency1 = []
for i in range(len(frequency1)):
    relative_frequency1.append(frequency1[i] / 200)
sumfrequency1 = sum(frequency1)
sumrelative_frequency1 = sum(relative_frequency1)

pan = pd.Series(data_norm)
frequency2 = pan.groupby(pd.cut(pan, bins=a2, right=True)).count()
frequency2 = frequency2.tolist()
frequency2[0] += 1
relative_frequency2 = []
for i in range(len(frequency2)):
    relative_frequency2.append(frequency2[i] / 200)
sumfrequency2 = sum(frequency2)
sumrelative_frequency2 = sum(relative_frequency2)

pan = pd.Series(data_unif)
frequency3 = pan.groupby(pd.cut(pan, bins=a3, right=True)).count()
frequency3 = frequency3.tolist()
relative_frequency3 = []
for i in range(len(frequency3)):
    relative_frequency3.append(frequency3[i] / 200)
sumfrequency3 = sum(frequency3)
sumrelative_frequency3 = sum(relative_frequency3)
# Задание 1
s_m1 = 0.0
for i in range(len(x1)):
    s_m1 += relative_frequency1[i] * x1[i]
labda1 = 1 / s_m1

f1 = []
F1 = []
for i in range(len(a1)):
    f1.append(labda1 * math.e ** (-labda1 * a1[i]))
    F1.append(1 - math.e ** (-labda1 * a1[i]))

p1 = []
for i in range(1, m):
    p1.append(F1[i] - F1[i - 1])
p1.append(1 - F1[7])
sump1 = sum(p1)

rz_fr1 = []
for i in range(len(relative_frequency1)):
    rz_fr1.append(abs(relative_frequency1[i] - p1[i]))

xi1 = []
for i in range(len(relative_frequency1)):
    xi1.append((200 * (relative_frequency1[i] - p1[i]) ** 2) / p1[i])
xi1sum = sum(xi1)
xi1crit = chi2.ppf(1 - .05, m - 2)

data = data_exp
h = (max(data_exp) - min(data_exp)) / m
x = np.arange(0, 3.5, 0.001)
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.hist(data, edgecolor='black', weights=(np.ones_like(data) / (len(data))) / h, bins=m)
y = scipy.stats.expon.pdf(x=x, scale=1/labda1)
plt.plot(x, y)
plt.savefig("myimage1.png", dpi=2000)
print(xi1sum, xi1crit)


# Задание 2
s_m1 = s_m2 = 0.0
for i in range(len(x2)):
    s_m1 += relative_frequency2[i] * x2[i]
    s_m2 += relative_frequency2[i] * x2[i] * x2[i]
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)
labda2 = 1 / s_m1

f2 = []
F2 = []
for i in range(len(a2)):
    f2.append(labda2 * math.e ** (-labda2 * a2[i]))
    F2.append(1 - math.e ** (-labda2 * a2[i]))
t1 = []
t2 = []
t3 = []
for i in range(len(a2)):
    t1.append((a2[i] - s_m1) / Smqd)
    t2.append((1 / Smqd) * (1 / math.sqrt(2 * math.pi) * math.e ** (-(t1[i] * t1[i] / 2))))
    t3.append(norm.cdf(t1[i]))

p2 = [t3[1]]
for i in range(2, m):
    p2.append(t3[i] - t3[i - 1])
p2.append(1 - t3[7])
print((1 / Smqd) * (1 / math.sqrt(2 * math.pi) * math.e ** (-(-0.08805  * -0.08805 / 2))))
sump2 = sum(p2)

rz_fr2 = []
for i in range(len(relative_frequency2)):
    rz_fr2.append(abs(relative_frequency2[i] - p2[i]))

xi2 = []
for i in range(len(relative_frequency2)):
    xi2.append((200 * (relative_frequency2[i] - p2[i]) ** 2) / p2[i])
xi2sum = sum(xi2)
xi2crit = chi2.ppf(1 - .05, m - 3)

data = data_norm
h = (max(data_norm) - min(data_norm)) / m
x = np.arange(-2, 5, 0.001)
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.hist(data, edgecolor='black', weights=(np.ones_like(data) / (len(data))) / h, bins=m)
plt.plot(x, norm.pdf(x, s_m1, Smqd))
plt.savefig("myimage2.png", dpi=2000)
print(xi2sum, xi2crit)
print(s_m1, Sample_variance)

# Задание 3
f3 = []
F3 = []
for i in range(len(a3)):
    f3.append(1 / (b - a))
    F3.append((a3[i] - a) / (b - a))
p3 = []
for i in range(1, m + 1):
    p3.append(F3[i] - F3[i - 1])
sump3 = sum(p3)

rz_fr3 = []
for i in range(len(relative_frequency3)):
    rz_fr3.append(abs(relative_frequency3[i] - p3[i]))

xi3 = []
for i in range(len(relative_frequency3)):
    xi3.append((200 * (relative_frequency3[i] - p3[i]) ** 2) / p3[i])
xi3sum = sum(xi3)
xi3crit = chi2.ppf(1 - .05, m - 1)

data = data_unif
x = np.arange(a, b, 0.001)
h = (max(data_unif) - min(data_unif)) / m
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
ax.hist(data, edgecolor='black', weights=(np.ones_like(data) / (len(data))) / h, bins=m)
y = scipy.stats.uniform.pdf(x=x, loc=a, scale=b - a)
plt.plot(x, y)
plt.savefig("myimage3.png", dpi=2000)
print(xi3sum, xi3crit)

# Задание 4
Dn4 = 0
y1 = []
for i in range(len(data_unif)):
    y1.append((data_unif[i] - a) / (b - a))
    if max(abs((i+1)/ 200 - (data_unif[i] - a) / (b - a)), abs((i ) / 200 - (data_unif[i] - a) / (b - a))) > Dn4:
        Dn4 = max(abs((i+1) / 200 - (data_unif[i] - a) / (b - a)), abs((i ) / 200 - (data_unif[i] - a) / (b - a)))
        xS4 = data_unif[i]
        j4 = i+1
if Dn4 * math.sqrt(200) < 1 - 0.05:
    test4 = True
else:
    test4 = False
answer4 = [a, b, 200, Dn4, Dn4 * math.sqrt(200), xS4, (xS4 - a) / (b - a), j4 / 200, (j4 - 1) / 200, test4]
x = data_unif
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
y = np.arange(1/200, 201/200, 1/200)
plt.plot(x, y)
plt.plot(x, y1)
plt.savefig("myimage4.png", dpi=2000)

# Задание 5
lambd = 1.53
Dn5 = 0
y2 =[]
for i in range(len(data_exp)):
    y2.append((1 - math.e ** (-lambd * data_exp[i])))
    if max(abs((i+1) / 200 - (1 - math.e ** (-lambd * data_exp[i]))), abs((i ) / 200 - (1 - math.e ** (-lambd * data_exp[i])))) > Dn5:
        Dn5 = max(abs((i+1) / 200 - (1 - math.e ** (-lambd * data_exp[i]))),
                  abs((i ) / 200 - (1 - math.e ** (-lambd * data_exp[i]))))
        xS5 = data_exp[i]
        j5 = i+1
if Dn5 * math.sqrt(200) < 1 - 0.05:
    test5 = True
else:
    test5 = False
answer5 = [lambd, 200, Dn5, Dn5 * math.sqrt(200), xS5, (1 - math.e ** (-lambd * data_exp[j5-1])),j5 / 200, (j5 - 1) / 200, test5]
print()
print(answer5)
print(answer4)
x = data_exp
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
y = np.arange(1/200, 201/200, 1/200)
plt.plot(x, y)
plt.plot(x, y2)
plt.show()
plt.savefig("myimage5.png", dpi=2000)

# Задание 1
f.write(str("Задание 1"))
f.write('\n')
f.write(str([round(elem, 5) for elem in data_exp]))
f.write('\n')
f.write(str([round(elem, 5) for elem in a1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in frequency1]))
f.write('\n')
f.write(str(round(sumfrequency1, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency1]))
f.write('\n')
f.write(str(round(sumrelative_frequency1, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in f1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in F1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in p1]))
f.write('\n')
f.write(str(round(sump1, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency1]))
f.write('\n')
f.write(str(round(sumrelative_frequency1, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in p1]))
f.write('\n')
f.write(str(round(sump1, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in rz_fr1]))
f.write('\n')
f.write(str(round(max(rz_fr1), 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in xi1]))
f.write('\n')
f.write(str(round(xi1sum, 5)))
f.write('\n')
f.write('\n')

# Задание 2
f.write(str("Задание 2"))
f.write('\n')
f.write(str([round(elem, 5) for elem in data_norm]))
f.write('\n')
f.write(str([round(elem, 5) for elem in a2]))
f.write('\n')
f.write(str([round(elem, 5) for elem in frequency2]))
f.write('\n')
f.write(str(round(sumfrequency2, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency2]))
f.write('\n')
f.write(str(round(sumrelative_frequency2, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a2]))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')
f.write(str([round(elem, 5) for elem in t2]))
f.write('\n')
f.write(str([round(elem, 5) for elem in t3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in p2]))
f.write('\n')
f.write(str(round(sump2, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a2]))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency2]))
f.write('\n')
f.write(str(round(sumrelative_frequency2, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in p2]))
f.write('\n')
f.write(str(round(sump2, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in rz_fr2]))
f.write('\n')
f.write(str(round(max(rz_fr2), 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in xi2]))
f.write('\n')
f.write(str(round(xi2sum, 5)))
f.write('\n')
f.write('\n')

# Задание 3
f.write(str("Задание 3"))
f.write('\n')
f.write(str([round(elem, 5) for elem in data_unif]))
f.write('\n')
f.write(str([round(elem, 5) for elem in a3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in frequency3]))
f.write('\n')
f.write(str(round(sumfrequency3, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency3]))
f.write('\n')
f.write(str(round(sumrelative_frequency3, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in f3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in F3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in p3]))
f.write('\n')
f.write(str(round(sump3, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in a3]))
f.write('\n')
f.write(str([round(elem, 5) for elem in relative_frequency3]))
f.write('\n')
f.write(str(round(sumrelative_frequency3, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in p3]))
f.write('\n')
f.write(str(round(sump3, 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in rz_fr3]))
f.write('\n')
f.write(str(round(max(rz_fr3), 5)))
f.write('\n')
f.write(str([round(elem, 5) for elem in xi3]))
f.write('\n')
f.write(str(round(xi3sum, 5)))
f.write('\n')
f.write('\n')

# Задание 4
f.write(str("Задание 4"))
f.write('\n')
f.write(str([round(elem, 5) for elem in answer4]))
f.write('\n')
f.write('\n')

# Задание 5
f.write(str("Задание 5"))
f.write('\n')
f.write(str([round(elem, 5) for elem in answer5]))
f.write('\n')
f.write('\n')

print(scipy.special.kolmogi(0.05))