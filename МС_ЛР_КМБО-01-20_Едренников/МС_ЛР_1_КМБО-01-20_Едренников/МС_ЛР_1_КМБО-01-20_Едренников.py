import numpy as np
import math
from scipy.stats import hypergeom
from scipy.stats import binom
from scipy.stats import norm
from scipy.stats import expon
from scipy.stats import uniform
import statistics
import pandas as pd

f = open('answer.txt', 'r+')
selection = np.random.binomial(n=11, p=0.23, size=200)
n = 11
p = 0.23
q = 0.77
selection_list = list(selection)
selection_numbers = list(set(selection))
frequency = []
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for number in selection_numbers:
    frequency.append(selection_list.count(number))
    relative_frequency.append(selection_list.count(number) / 200)
    ter_frequency.append(binom.pmf(k=number, n=11, p=0.23))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

sum_frequency = []
for _ in range(len(relative_frequency)):
    sum_frequency.append(None)

for i in range(len(relative_frequency)):
    sum_frequency[i] = 0
    for j in range(i):
        sum_frequency[i] += relative_frequency[j]

mean = 0.0
for i in range(len(selection_numbers)):
    mean += frequency[i] * selection_numbers[i]
mean = mean / 200

s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(selection_numbers)):
    s_m1 += relative_frequency[i] * selection_numbers[i]
    s_m2 += relative_frequency[i] * (selection_numbers[i] ** 2)
    s_m3 += relative_frequency[i] * (selection_numbers[i] ** 3)
    s_m4 += relative_frequency[i] * (selection_numbers[i] ** 4)
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)

mode = statistics.multimode(selection_list)

median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], mode[0], median, S_a_c, S_k_c]
last_table_T = [n * p, n * p * q, math.sqrt(n * p * q), math.floor((n + 1) * p), (n + 1) * p - 0.5, math.ceil(n * p),
                (q - p) / math.sqrt(n * p * q), (1 - 6 * p * q) / (n * p * q)]
L_rz = []
for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)

f.write(str(selection))
f.write('\n')
selection_list.sort()
f.write(str(selection_list))
f.write('\n')
f.write(str(selection_numbers))
f.write('\n')
f.write(str(frequency))
f.write('\n')
f.write(str(relative_frequency))
f.write('\n')
f.write(str(sum_frequency))
f.write('\n')
f.write(str(ter_frequency))
f.write('\n')
f.write(str(rz))
f.write('\n')
f.write(str(last_table))
f.write('\n')
f.write(str(last_table_T))
f.write('\n')
f.write(str(L_rz))
f.write('\n')
f.write(str(O_rz))
f.write('\n')
f.write('\n')

selection = np.random.geometric(p=0.23, size=200)
selection = [0,	1,	0,	0,	2,	0,	9,	6,	0,	18,
0,	0,	0,	2,	5,	11,	7,	5,	2,	4,
5,	7,	2,	1,	0,	0,	0,	0,	0,	2,
5,	5,	1,	1,	5,	2,	19,	1,	1,	1,
14,	0,	1,	1,	3,	0,	2,	1,	0,	10,
7,	0,	15,	0,	1,	7,	3,	2,	1,	0,
0,	3,	4,	4,	3,	6,	2,	1,	2,	2,
0,	3,	5,	1,	9,	1,	2,	1,	2,	1,
1,	0,	6,	6,	2,	0,	4,	1,	0,	6,
5,	3,	4,	1,	20,	0,	2,	8,	3,	2,
3,	1,	0,	0,	10,	0,	3,	0,	5,	2,
7,	1,	3,	0,	6,	6,	7,	0,	2,	0,
7,	4,	2,	6,	0,	2,	8,	3,	0,	0,
7,	2,	10,	1,	0,	11,	4,	0,	1,	2,
6,	0,	3,	2,	3,	6,	0,	3,	0,	0,
1,	0,	0,	3,	9,	4,	8,	1,	0,	4,
4,	13,	1,	6,	15,	1,	10,	0,	1,	0,
0,	0,	2,	4,	1,	11,	1,	10,	0,	13,
1,	4,	6,	1,	1,	0,	6,	0,	3,	1,
3,	6,	8,	1,	4,	0,	1,	1,	5,	1,]
q = 0.77
p = 0.23
selection_list = list(selection)
selection_numbers = list(set(selection))
frequency = []
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for number in selection_numbers:
    frequency.append(selection_list.count(number))
    relative_frequency.append(selection_list.count(number) / 200)
    ter_frequency.append((q ** (number - 1)) * p)

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

sum_frequency = []
for _ in range(len(relative_frequency)):
    sum_frequency.append(None)

for i in range(len(relative_frequency)):
    sum_frequency[i] = 0
    for j in range(i):
        sum_frequency[i] += relative_frequency[j]

mean = 0.0
for i in range(len(selection_numbers)):
    mean += frequency[i] * selection_numbers[i]
mean = mean / 200

s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(selection_numbers)):
    s_m1 += relative_frequency[i] * selection_numbers[i]
    s_m2 += relative_frequency[i] * (selection_numbers[i] ** 2)
    s_m3 += relative_frequency[i] * (selection_numbers[i] ** 3)
    s_m4 += relative_frequency[i] * (selection_numbers[i] ** 4)
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)

mode = statistics.multimode(selection_list)

median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, median, S_a_c, S_k_c]
last_table_T = [q / p, q / (p * p), math.sqrt(q) / p, 0, math.floor(-math.log(2) / math.log(q)),
                -math.log(2) / math.log(q) - 0.5, (2 - p) / math.sqrt(q), 6 + (p * p) / q]
L_rz = []

for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)

f.write(str(selection))
f.write('\n')
selection_list.sort()
f.write(str(selection_list))
f.write('\n')
f.write(str(selection_numbers))
f.write('\n')
f.write(str(frequency))
f.write('\n')
f.write(str(relative_frequency))
f.write('\n')
f.write(str(sum_frequency))
f.write('\n')
f.write(str(ter_frequency))
f.write('\n')
f.write(str(rz))
f.write('\n')
f.write(str(last_table))
f.write('\n')
f.write(str(last_table_T))
f.write('\n')
f.write(str(L_rz))
f.write('\n')
f.write(str(O_rz))
f.write('\n')
f.write('\n')

selection = np.random.poisson(0.56, 200)
p = 0.56
selection_list = list(selection)
selection_numbers = list(set(selection))
frequency = []
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for number in selection_numbers:
    frequency.append(selection_list.count(number))
    relative_frequency.append(selection_list.count(number) / 200)
    ter_frequency.append(((p ** number) / (math.factorial(number))) * math.e ** (-p))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

sum_frequency = []
for _ in range(len(relative_frequency)):
    sum_frequency.append(None)

for i in range(len(relative_frequency)):
    sum_frequency[i] = 0
    for j in range(i):
        sum_frequency[i] += relative_frequency[j]

mean = 0.0
for i in range(len(selection_numbers)):
    mean += frequency[i] * selection_numbers[i]
mean = mean / 200

s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(selection_numbers)):
    s_m1 += relative_frequency[i] * selection_numbers[i]
    s_m2 += relative_frequency[i] * (selection_numbers[i] ** 2)
    s_m3 += relative_frequency[i] * (selection_numbers[i] ** 3)
    s_m4 += relative_frequency[i] * (selection_numbers[i] ** 4)
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)

mode = statistics.multimode(selection_list)

median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, S_a_c, S_k_c]
last_table_T = [p, p, math.sqrt(p), math.floor(p), math.floor(p + 1 / 3 - 0.02 / p), 1 / math.sqrt(p),
                1 / p]
L_rz = []

for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)

f.write(str(selection))
f.write('\n')
selection_list.sort()
f.write(str(selection_list))
f.write('\n')
f.write(str(selection_numbers))
f.write('\n')
f.write(str(frequency))
f.write('\n')
f.write(str(relative_frequency))
f.write('\n')
f.write(str(sum_frequency))
f.write('\n')
f.write(str(ter_frequency))
f.write('\n')
f.write(str(rz))
f.write('\n')
f.write(str(last_table))
f.write('\n')
f.write(str(last_table_T))
f.write('\n')
f.write(str(L_rz))
f.write('\n')
f.write(str(O_rz))
f.write('\n')
f.write('\n')
selection = np.random.random_integers(0, 10, 200)
n = 11
selection_list = list(selection)
selection_numbers = list(set(selection))
frequency = []
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for number in selection_numbers:
    frequency.append(selection_list.count(number))
    relative_frequency.append(selection_list.count(number) / 200)
    ter_frequency.append(1 / 11)

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

sum_frequency = []
for _ in range(len(relative_frequency)):
    sum_frequency.append(None)

for i in range(len(relative_frequency)):
    sum_frequency[i] = 0
    for j in range(i):
        sum_frequency[i] += relative_frequency[j]

mean = 0.0
for i in range(len(selection_numbers)):
    mean += frequency[i] * selection_numbers[i]
mean = mean / 200

s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(selection_numbers)):
    s_m1 += relative_frequency[i] * selection_numbers[i]
    s_m2 += relative_frequency[i] * (selection_numbers[i] ** 2)
    s_m3 += relative_frequency[i] * (selection_numbers[i] ** 3)
    s_m4 += relative_frequency[i] * (selection_numbers[i] ** 4)
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)

mode = statistics.multimode(selection_list)

median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, S_a_c, S_k_c]
last_table_T = [(n - 1) / 2, ((n * n) - 1) / 12, 0.5 * math.sqrt(((n * n) - 1) / 3), (n - 1) / 2, (n - 1) / 2,
                0, -(6 / 5) * (((n * n) + 1) / ((n * n)) - 1)]
L_rz = []

for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)

f.write(str(selection))
f.write('\n')
selection_list.sort()
f.write(str(selection_list))
f.write('\n')
f.write(str(selection_numbers))
f.write('\n')
f.write(str(frequency))
f.write('\n')
f.write(str(relative_frequency))
f.write('\n')
f.write(str(sum_frequency))
f.write('\n')
f.write(str(ter_frequency))
f.write('\n')
f.write(str(rz))
f.write('\n')
f.write(str(last_table))
f.write('\n')
f.write(str(last_table_T))
f.write('\n')
f.write(str(L_rz))
f.write('\n')
f.write(str(O_rz))
f.write('\n')
f.write('\n')
selection = np.random.hypergeometric(12, 13, 11, 200)
m = 11
M = 25
K = 12
selection_list = list(selection)
selection_numbers = list(set(selection))
frequency = []
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for number in selection_numbers:
    frequency.append(selection_list.count(number))
    relative_frequency.append(selection_list.count(number) / 200)
    ter_frequency.append(hypergeom.cdf(number, 25, 11, 12))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

sum_frequency = []
for _ in range(len(relative_frequency)):
    sum_frequency.append(None)

for i in range(len(relative_frequency)):
    sum_frequency[i] = 0
    for j in range(i):
        sum_frequency[i] += relative_frequency[j]

mean = 0.0
for i in range(len(selection_numbers)):
    mean += frequency[i] * selection_numbers[i]
mean = mean / 200

s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(selection_numbers)):
    s_m1 += relative_frequency[i] * selection_numbers[i]
    s_m2 += relative_frequency[i] * (selection_numbers[i] ** 2)
    s_m3 += relative_frequency[i] * (selection_numbers[i] ** 3)
    s_m4 += relative_frequency[i] * (selection_numbers[i] ** 4)
Sample_variance = s_m2 - (s_m1 ** 2)
Smqd = math.sqrt(Sample_variance)

mode = statistics.multimode(selection_list)

median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, median, S_a_c, S_k_c]
last_table_T = [m * K / M, (m * K * (M - K) * (M - m)) / ((M - 1) * (M * M)),
                (1 / M) * math.sqrt((m * K * (M - K) * (M - m)) / (M - 1)),
                math.floor(((K + 1) * (m - 1)) / (M + 2)), m, m + 0.5,
                ((M - 2 * K) * (M - 2 * m) / (M - 2)) * math.sqrt((M - 1) / (m * K * (M - K) * (M - m))),
                math.floor(((M - 1) * M * M) / (m * (M - 2) * (M - 3) * (M - m))) * math.floor(
                    (M * (M + 1) - 6 * M * (M - m)) / (K * (M - K)) + (3 * m * (M + 6) * (M - m)) / (M * M) - 6)]
L_rz = []

for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)

f.write(str(selection))
f.write('\n')
selection_list.sort()
f.write(str(selection_list))
f.write('\n')
f.write(str(selection_numbers))
f.write('\n')
f.write(str(frequency))
f.write('\n')
f.write(str(relative_frequency))
f.write('\n')
f.write(str(sum_frequency))
f.write('\n')
f.write(str(ter_frequency))
f.write('\n')
f.write(str(rz))
f.write('\n')
f.write(str(last_table))
f.write('\n')
f.write(str(last_table_T))
f.write('\n')
f.write(str(L_rz))
f.write('\n')
f.write(str(O_rz))
f.write('\n')
f.write('\n')
f.close()

f = open('answer1.txt', 'r+')
selection = np.random.normal(0.3, 1.03 , 200)
selection = [1.1658, 0.03992, 1.01534, 0.07553, -0.92043, 0.70086, 0.35958, 0.39276, -1.48112, -1.2337, -1.30807, -2.14042, 0.45475, -1.31484, -0.4506, -0.69883, 1.72924, -0.71028, -0.39715, 0.81071, -0.6385, 0.40144, 0.92914, 0.6065, 1.07279, 0.74965, 1.60671, 2.28143, -0.07828, -0.1886, 0.64163, -0.45412, 0.89457, 4.23547, 1.1833, -0.13018, -0.51238, -0.61987, 0.62851, 0.23584, 3.21512, -1.08879, 1.81319, -0.34664, -0.05514, -1.39536, 0.85688, 0.71801, -0.31807, 2.36985, -0.43662, 0.18998, -0.17463, 0.46352, 0.60267, 1.0832, -0.55948, 0.15433, -1.33625, -1.15913, 1.23691, -1.02543, -0.47616, 0.99325, 1.10945, -1.32368, 0.31241, 0.45591, 0.65407, -0.53396, 1.04851, -2.02812, -2.1909, 1.68378, 0.33701, 0.15507, -0.64901, 0.23099, -1.00966, 0.80459, -1.05305, -0.13572, 0.29378, 0.01389, 0.27026, 0.82507, 1.00571, 2.0791, 0.16196, 1.2796, 0.86837, 0.11846, -0.36415, -0.20819, 0.14569, 0.67373, -1.22727, -0.8153, -0.03751, 0.92959, 0.01598, 0.13447, -0.82175, 0.95921, 1.02299, 0.89841, -1.11187, 0.47599, -0.117, -0.73716, -0.80665, -0.83037, -0.90396, -1.93622, 0.88252, 0.87349, 0.84677, 0.19989, 0.17481, -0.5227, 0.59634, 0.56526, 0.75848, -0.31772, -1.12639, 0.40775, 0.6983, 0.95638, 1.16627, 1.30625, 1.00296, 0.90736, 0.72462, 1.44972, 0.86584, -1.35144, -1.11026, 0.39271, 1.04747, -1.65366, 0.37406, -0.65462, -0.0471, 0.62038, 1.40446, 0.97301, -0.62792, -1.89807, 0.04357, 0.34359, -0.03187, 0.91351, 0.8047, 1.29796, -0.65203, 0.44832, 0.87969, -1.10788, 1.38339, -0.86818, -0.05037, -1.73752, -0.0647, -0.43745, 0.3148, -0.52213, 1.62558, 0.2308, 2.37026, 0.40379, 0.64067, -0.49419, 0.3358, 0.85429, -0.61423, 1.34232, 0.94634, 0.33478, -0.75317, 0.82504, -0.91629, -1.78466, -0.6901, 0.98902, 1.65466, 2.0217, -2.10385, -1.0686, -0.51435, -0.24632, 1.07106, -0.2019, 1.37149, 1.14914, -0.04509, 0.02509, 0.18447, 1.84519, 0.30137, 0.24753]


mu = 0.3
sigma = 1.03 ** 2
selection_list = list(selection)
m = 1 + math.floor(math.log2(200))
a = [None] * (m+1)
x = []
a[0] = min(selection_list)
a[m] = max(selection_list)
for i in range(1, m):
    a[i] = a[i - 1] + (a[m] - a[0]) / m
a.sort()
for i in range(1, m):
    x.append((a[i - 1] + a[i]) / 2)


pan = pd.Series(selection_list)
print(pan)
frequency = pan.groupby(pd.cut(pan, bins=a, right=True)).count()
frequency = frequency.tolist()
frequency[0]+=1
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for i in range(len(frequency)):
    relative_frequency.append(frequency[i] / 200)
    ter_frequency.append(abs(norm.cdf(a[i], loc=mu, scale=1.03) - norm.cdf(a[i + 1], loc=mu, scale=1.03)))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

mean = 0.0
for i in range(len(x)):
    mean += frequency[i] * x[i]
mean = mean / 200
s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(x)):
    s_m1 += relative_frequency[i] * x[i]
    s_m2 += relative_frequency[i] * (x[i] ** 2)
    s_m3 += relative_frequency[i] * (x[i] ** 3)
    s_m4 += relative_frequency[i] * (x[i] ** 4)
Sample_variance = 0.0
h = (a[m-1]-a[0])/m
for i in range(1, m-1):
    Sample_variance += ((x[i] - mean)**2)*relative_frequency[i]
Sample_variance = Sample_variance - (h*h)/12
Smqd = math.sqrt(Sample_variance)

k = relative_frequency.index(max(relative_frequency))
ak = a[k]
mode = [ak + h * (relative_frequency[k] - relative_frequency[k - 1]) / (
2 * relative_frequency[k] - relative_frequency[k - 1] - relative_frequency[
(k + 1) % len(relative_frequency)])]
median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, S_a_c, S_k_c]
last_table_T = [0.3, 1.03 ** 2, 1.03, 0.3, 0.3, 0, 0]
L_rz = []
for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)
f.write(str([ round(elem, 5) for elem in selection ]))
f.write('\n')
selection_list.sort()
f.write(str([ round(elem, 5) for elem in selection_list ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in a ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in x ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in relative_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in ter_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table_T ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in L_rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in O_rz ]))
f.write('\n')
f.write('\n')


selection = expon.rvs(scale=1/2.06, size=200)
selection = [0.37721, 0.03169, 0.37151, 0.0769, 0.34037, 1.13088, 0.2575, 0.56454, 1.82127, 0.05064, 0.10365, 0.70715, 0.21093, 0.37443, 0.19333, 0.09039, 0.04419, 0.12236, 0.77516, 0.40823, 1.4607, 1.83618, 1.37025, 0.25452, 0.68557, 1.03029, 1.14552, 0.2294, 0.28135, 0.73949, 0.07434, 0.10757, 1.42173, 0.07836, 0.04259, 0.01033, 0.21435, 0.38462, 0.9403, 1.08397, 0.95688, 0.10289, 0.0546, 1.92442, 1.40552, 0.09963, 0.59998, 0.12916, 0.11015, 0.39198, 0.52947, 0.03863, 0.1736, 0.10894, 1.18496, 0.12376, 1.5779, 0.68347, 0.66878, 1.07169, 1.08322, 0.4471, 1.59374, 0.04513, 0.16902, 0.67763, 0.25521, 0.359, 0.34856, 0.21629, 2.68677, 1.53046, 0.24074, 0.15679, 0.18262, 0.12134, 1.59947, 1.40134, 0.83365, 1.5909, 0.00551, 0.16201, 0.19019, 0.46259, 1.81955, 0.03282, 0.98149, 0.56862, 0.49926, 0.65135, 0.0878, 0.06003, 1.27631, 0.02527, 0.12644, 0.53766, 0.44759, 0.01381, 0.75243, 0.096, 0.43248, 0.079, 0.06444, 0.08693, 0.05088, 0.62982, 0.66711, 0.00541, 0.28463, 0.1487, 0.63705, 0.19284, 0.03226, 0.38898, 1.77783, 0.3157, 0.05891, 0.54325, 0.71747, 1.06411, 0.32912, 0.34705, 0.55192, 0.01571, 0.44481, 0.43602, 0.59003, 0.24869, 0.1213, 0.16289, 0.46451, 0.59355, 0.35801, 1.29891, 0.18403, 0.15223, 2.21691, 1.21771, 0.02043, 0.00464, 0.17036, 1.13631, 0.33452, 0.09827, 2.42315, 0.94445, 0.25844, 0.36272, 0.56341, 0.19087, 0.1771, 0.26055, 0.09611, 0.48435, 0.09605, 0.57185, 0.46519, 1.0077, 0.18623, 0.2779, 0.08348, 0.01262, 0.05042, 0.47721, 0.36026, 0.74537, 0.05519, 0.49033, 1.72484, 0.15356, 0.30864, 0.19948, 1.60141, 0.00222, 0.0428, 1.12284, 0.75679, 0.30892, 0.00987, 0.91402, 0.01588, 0.1943, 0.14746, 0.35594, 1.59502, 1.02987, 0.4608, 0.60957, 0.19554, 0.74813, 0.01426, 0.11249, 1.66962, 0.44219, 0.71825, 1.16398, 0.04346, 0.45743, 0.60865, 0.14083]

lambd = 2.06
selection_list = list(selection)
m = 1 + math.floor(math.log2(200))
a = [None] * (m+1)
x = []
a[0] = min(selection_list)
a[m] = max(selection_list)

a[0] = 0
a[m] = max(selection_list)
for i in range(1, m):
    a[i] = a[i - 1] + (a[m] - a[0]) / m
for i in range(1, m):
    x.append((a[i - 1] + a[i]) / 2)

pan = pd.Series(selection_list)
frequency = pan.groupby(pd.cut(pan, bins=a, right=True)).count()
frequency = frequency.tolist()
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []
print(sum(frequency))
for i in range(len(frequency)):
    relative_frequency.append(frequency[i] / 200)
    ter_frequency.append(abs(expon.cdf(x = a[i], scale = 1/lambd) - expon.cdf(x = a[i+1], scale = 1/lambd)))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

mean = 0.0
for i in range(len(x)):
    mean += frequency[i] * x[i]
mean = mean / 200
s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(x)):
    s_m1 += relative_frequency[i] * x[i]
    s_m2 += relative_frequency[i] * (x[i] ** 2)
    s_m3 += relative_frequency[i] * (x[i] ** 3)
    s_m4 += relative_frequency[i] * (x[i] ** 4)
Sample_variance = 0.0
h = (a[m-1]-a[0])/m
for i in range(1, m-1):
    Sample_variance += ((x[i] - mean)**2)*relative_frequency[i]
Sample_variance = Sample_variance - (h*h)/12
Smqd = math.sqrt(Sample_variance)

k = relative_frequency.index(max(relative_frequency))
ak = a[k]
mode = [ak+h*(relative_frequency[k]-relative_frequency[len(relative_frequency)-1])/(2*relative_frequency[k]-relative_frequency[len(relative_frequency)-1]-relative_frequency[(k+1)%len(relative_frequency)])]
median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, S_a_c, S_k_c]
last_table_T = [lambd**(-1), lambd**(-2), lambd**(-1), 0, math.log(2)/lambd, 2, 6]
L_rz = []

for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)
f.write(str([ round(elem, 5) for elem in selection ]))
f.write('\n')
selection_list.sort()
f.write(str([ round(elem, 5) for elem in selection_list ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in a ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in x ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in relative_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in ter_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table_T ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in L_rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in O_rz ]))
f.write('\n')
f.write('\n')

selection = np.random.uniform(0.12, 6.12, 200)
selection = [0.12562551492260876, 0.1272326588196918, 0.1369396295741897, 0.17533268565879212, 0.1807999305691469, 0.1921757167092718, 0.19922495406104213, 0.24435789958973364, 0.34443804186878835, 0.3995417773692066, 0.42734823985142023, 0.4329020890151357, 0.46189051602031606, 0.48988127406484605, 0.5072492218424186, 0.628531535239479, 0.6361746981421076, 0.6479075503856279, 0.6902630579657861, 0.7076354415522476, 0.7425551559881848, 0.7645196884125206, 0.7931707333494434, 0.7987184351648705, 0.811048132234817, 0.8507106637667311, 0.86442346827153, 0.876692275639506, 0.9305788546920793, 0.9546073496588595, 0.9781384063000343, 0.9901717653170478, 0.992767961783935, 0.9950895358772692, 1.0278481166156426, 1.052054688287091, 1.1182467556776534, 1.153118323529029, 1.2054330726271774, 1.2484579728911331, 1.282558646818837, 1.3008212154995515, 1.3159847155442832, 1.3436283386544376, 1.3595550986589968, 1.3602380389897575, 1.3750705715803093, 1.3975777101672864, 1.398207992944255, 1.4003853099675103, 1.4246849631058036, 1.5219092875121407, 1.5275875048183982, 1.5952994081713032, 1.6350676382701672, 1.6400779088455586, 1.68066908371132, 1.701360573842325, 1.7068232222260717, 1.7181919302738398, 1.7506514707796978, 1.7530338313747982, 1.7643062656788473, 1.7832642822541285, 1.839229338785609, 1.8986498908953924, 1.9327810326298858, 1.9898800135065988, 1.9979161756900243, 2.0017580685950485, 2.046079510022451, 2.1516443376470162, 2.2269025853612874, 2.2531484934685, 2.2658406730212564, 2.2908572477682, 2.3891815415675968, 2.401152044373532, 2.405267589004157, 2.4158320858197246, 2.422384286776534, 2.556404467799983, 2.565914416211781, 2.595107396043554, 2.6358160002948114, 2.652512980596657, 2.658273841307581, 2.668717784784176, 2.677177824325412, 2.6777581700601427, 2.6928353019824725, 2.776512530004556, 2.807953558606706, 2.8106224456249422, 2.8674268022000646, 2.931440242786512, 3.029814884530647, 3.0574862508707166, 3.06776575268203, 3.107825604136611, 3.1712246494681198, 3.183444247937006, 3.19558988827834, 3.279580535680761, 3.34874262456642, 3.3653266087751876, 3.4130117547429895, 3.4289558467876464, 3.429707884597394, 3.4612265217230136, 3.472472681463942, 3.4931555108173145, 3.520725244369217, 3.538771657082056, 3.5419189052203004, 3.5598404007997533, 3.56077259680631, 3.6025890658399033, 3.6184842317147483, 3.622598963028352, 3.6312163624616987, 3.673844689419713, 3.699360311574848, 3.699584238184599, 3.7394901362126474, 3.743093645304036, 3.7466573687129108, 3.7528664752726737, 3.759502542199381, 3.7635750497138556, 3.7759918444216707, 3.8190581090850637, 3.8430687878151777, 3.8649772373242564, 3.8927232051271354, 3.924090156694729, 3.9259508707491255, 3.929232585775269, 3.9858348504115755, 3.9889380803560135, 4.008942857938251, 4.025138255026329, 4.040541262518161, 4.0923095856343705, 4.095987037332115, 4.098072389259218, 4.143618301255864, 4.155121974493067, 4.278022360228854, 4.29167788368445, 4.317557031823923, 4.396709002172995, 4.428118100839357, 4.432551318956099, 4.465796351565596, 4.473200322801036, 4.491543410966609, 4.516186991997199, 4.565216838242331, 4.59439355855376, 4.603669341946267, 4.6595402461291275, 4.663275676373524, 4.687657717534638, 4.711675845002831, 4.743835306041133, 4.793160932614812, 4.796558074792212, 4.9125812290668565, 4.922395516845044, 4.981580640571454, 5.106257268972576, 5.132195235230196, 5.1888707752868, 5.29894406056668, 5.329265876480405, 5.401803167875951, 5.4326323551708935, 5.512847023208698, 5.540167462105237, 5.573194892165703, 5.633699796719307, 5.657349195565742, 5.658525137804989, 5.740116246629289, 5.74956982419598, 5.80176546829411, 5.8120331764405435, 5.839235088038437, 5.849153649751552, 5.871716631227106, 5.933824017837832, 5.9375642131599005, 5.944084152326927, 5.973435517157857, 5.983455810514716, 6.010887918991185, 6.046756999572009, 6.103750055130712, 6.106887898494262]

a1 = 0.12
b = 6.12
selection_list = list(selection)
m = 1 + math.floor(math.log2(200))
a = [None] * (m+1)
x = []
a[0] = a1
a[m] = b

for i in range(1, m ):
    a[i] = a[i - 1] + (a[m] - a[0]) / m
a.sort()
for i in range(1, m):
    x.append((a[i - 1] + a[i]) / 2)

pan = pd.Series(selection_list)
frequency = pan.groupby(pd.cut(pan, bins=a, right=True)).count()
frequency = frequency.tolist()
relative_frequency = []
ter_frequency = []
rz = []
O_rz = []

for i in range(len(frequency)):
    relative_frequency.append(frequency[i] / 200)
    ter_frequency.append(abs(uniform.cdf(x = a[i], loc = a1, scale = 6) - uniform.cdf(x = a[i+1], loc = a1, scale = 6)))

for i in range(len(relative_frequency)):
    rz.append(abs(relative_frequency[i] - ter_frequency[i]))

mean = 0.0
for i in range(len(x)):
    mean += frequency[i] * x[i]
mean = mean / 200
s_m1 = s_m2 = s_m3 = s_m4 = 0.0
for i in range(len(x)):
    s_m1 += relative_frequency[i] * x[i]
    s_m2 += relative_frequency[i] * (x[i] ** 2)
    s_m3 += relative_frequency[i] * (x[i] ** 3)
    s_m4 += relative_frequency[i] * (x[i] ** 4)
Sample_variance = 0.0
h = (a[m-1]-a[0])/m
for i in range(1, m-1):
    Sample_variance += ((x[i] - mean)**2)*relative_frequency[i]
Sample_variance = Sample_variance - (h*h)/12
Smqd = math.sqrt(Sample_variance)

k = relative_frequency.index(max(relative_frequency))
ak = a[k-1]
mode = [ak+h*(relative_frequency[k]-relative_frequency[k-1])/(2*relative_frequency[k]-relative_frequency[k-1]-relative_frequency[(k+1)%len(relative_frequency)])]
median = statistics.median(selection_list)

S_a_c = (s_m3 - 3 * s_m2 * s_m1 + 2 * (s_m1 ** 3)) / (Smqd ** 3)
S_k_c = (s_m4 - 4 * s_m3 * s_m1 + 6 * s_m2 * (s_m1 ** 2) - 3 * (s_m1 ** 4)) / (Smqd ** 4) - 3

last_table = [mean, Sample_variance, Smqd, mode[0], median, S_a_c, S_k_c]
last_table_T = [(a1+b)/2, ((a1+b)*(a1+b))/12, (b-a1)/(2*math.sqrt(3)), (a1+b)/2, (a1+b)/2, 0, -6/5]
L_rz = []
for i in range(len(last_table_T)):
    L_rz.append(abs(last_table[i] - last_table_T[i]))
    if last_table_T[i] != 0:
        O_rz.append(L_rz[i] / last_table_T[i])
    if last_table_T[i] == 0 and L_rz[i] == 0:
        O_rz.append(0)
f.write(str([ round(elem, 5) for elem in selection ]))
f.write('\n')
selection_list.sort()
f.write(str([ round(elem, 5) for elem in selection_list ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in a ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in x ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in relative_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in ter_frequency ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in last_table_T ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in L_rz ]))
f.write('\n')
f.write(str([ round(elem, 5) for elem in O_rz ]))
f.write('\n')
f.write('\n')
f.close()


