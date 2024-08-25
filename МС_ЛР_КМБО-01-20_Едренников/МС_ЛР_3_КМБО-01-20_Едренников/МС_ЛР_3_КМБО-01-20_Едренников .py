import math
import scipy.stats

f = open('answer.txt', 'r+')

"Задание 1"
x1 = [2.86728, 1.13188, 1.65331, -0.24485, 0.77905, 2.30288, 2.63247, -0.68905, 0.31365, 3.83315, 0.78376,
      0.53350, 2.64152, -0.41601, 2.47467, 0.87717, -0.04131, 1.26223, 0.79233, 1.58629]
x2 = [0.64330, -0.13516, 0.29311, -0.57104, 0.11235, 1.26531, 0.46257, 2.35009, 1.24173, 1.51367, 2.39267, 2.97226,
      1.18987, 0.49294, 1.34834, 2.37526, 0.53450, 1.61054, -1.60799, 1.50408]
x3 = [0.56871, 2.79170, -0.51887, 1.53098, 1.39604, -0.25179, 2.48092, 2.74350, -0.57746, 4.21158, 2.04421, 1.64747,
      3.30493, 3.57286, -0.24358, 1.65645, 0.36705, 1.76563, 4.39980, 2.46199]

x1s = 1 / 20 * sum(x1)
x2s = 1 / 20 * sum(x2)
x3s = 1 / 20 * sum(x3)

x12 = [x ** 2 for x in x1]
x22 = [x ** 2 for x in x2]
x32 = [x ** 2 for x in x3]

x1s2 = 1 / 20 * sum(x12)
x2s2 = 1 / 20 * sum(x22)
x3s2 = 1 / 20 * sum(x32)

S2x1 = (20 / 19) * (x1s2 - x1s)
S2x2 = (20 / 19) * (x2s2 - x2s)
S2x3 = (20 / 19) * (x3s2 - x3s)

T1 = ((x1s - x2s) / math.sqrt(19 * (S2x1 + S2x2))) * math.sqrt((20 * 20 * 38) / 40)
T2 = ((x1s - x3s) / math.sqrt(19 * (S2x1 + S2x3))) * math.sqrt((20 * 20 * 38) / 40)
T3 = ((x2s - x3s) / math.sqrt(19 * (S2x2 + S2x3))) * math.sqrt((20 * 20 * 38) / 40)

tcrit = scipy.stats.t.ppf(0.975, 38)

t11 = [x1s, x2s, x1s2, x2s2, S2x1, S2x2, T1]
t12 = [x1s, x3s, x1s2, x3s2, S2x1, S2x3, T2]
t13 = [x2s, x3s, x2s2, x3s2, S2x2, S2x3, T3]
t21 = []
t22 = []
t23 = []

if abs(T1) > tcrit:
    t21 = [abs(T1), tcrit, 0]
else:
    t21 = [abs(T1), tcrit, 1]

if abs(T2) > tcrit:
    t22 = [abs(T2), tcrit, 0]
else:
    t22 = [abs(T2), tcrit, 1]

if abs(T3) > tcrit:
    t23 = [abs(T3), tcrit, 0]
else:
    t23 = [abs(T3), tcrit, 1]

f.write(str("Задание 1"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t11]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t12]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t13]))
f.write('\n')

f.write(str("Таблица 2"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t21]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t22]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t23]))
f.write('\n')

"Задание 2"
xg = (1 / 3) * (x1s + x2s + x3s)
Sg = 0
for i in range(len(x1)):
    Sg += (x1[i] - xg) ** 2
    Sg += (x2[i] - xg) ** 2
    Sg += (x3[i] - xg) ** 2

Sf = 1 / 20 * ((x1s - xg) ** 2 + (x2s - xg) ** 2 + (x3s - xg) ** 2)
Sl = Sg - Sf

Sf2 = Sf / 2
Sl2 = Sl / 57

k1 = 2
k2 = 57

F = Sf2 / Sl2
Fcrit = scipy.stats.f.ppf(0.95, k1, k2)

t1 = [Sg, Sf, Sl, Sf2, Sl2, k1, k2, F]

t2 = []

if F > Fcrit:
    t2 = [F, 0.05, Fcrit, 0]
else:
    t2 = [F, 0.05, Fcrit, 1]

f.write(str("Задание 2"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')

f.write(str("Таблица 2"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t2]))
f.write('\n')

"Задание 3"
alpha = 0.05
t1, pval1 = scipy.stats.ttest_ind(x1, x2, equal_var=True)
t1, pval2 = scipy.stats.ttest_ind(x1, x3, equal_var=True)
t1, pval3 = scipy.stats.ttest_ind(x2, x3, equal_var=True)


if pval1 < alpha:
    t1 = [pval1, alpha, 0]
else:
    t1 = [pval1, alpha, 1]

if pval2 < alpha:
    t2 = [pval2, alpha, 0]
else:
    t2 = [pval2, alpha, 1]

if pval3 < alpha:
    t3 = [pval3, alpha, 0]
else:
    t3 = [pval3, alpha, 1]

f.write(str("Задание 3"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t2]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t3]))
f.write('\n')

"Задание 4"
alpha = 0.05
t1, pval1 = scipy.stats.ttest_ind(x1, x2, equal_var=False)
t1, pval2 = scipy.stats.ttest_ind(x1, x3, equal_var=False)
t1, pval3 = scipy.stats.ttest_ind(x2, x3, equal_var=False)

if pval1 < alpha:
    t1 = [pval1, alpha, 0]
else:
    t1 = [pval1, alpha, 1]

if pval2 < alpha:
    t2 = [pval2, alpha, 0]
else:
    t2 = [pval2, alpha, 1]

if pval3 < alpha:
    t3 = [pval3, alpha, 0]
else:
    t3 = [pval3, alpha, 1]

f.write(str("Задание 4"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t2]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t3]))
f.write('\n')

"Задание 5"
alpha = 0.05
t1, pval = scipy.stats.f_oneway(x1, x2, x3)
if pval < alpha:
    t1 = [pval, alpha, 0]
else:
    t1 = [pval, alpha, 1]

f.write(str("Задание 5"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')

"Задание 6"
Smax1 = max(S2x1, S2x2)
Smin1 = min(S2x1, S2x2)
Smax2 = max(S2x1, S2x3)
Smin2 = min(S2x1, S2x3)
Smax3 = max(S2x2, S2x3)
Smin3 = min(S2x2, S2x3)

k1 = k2 = 19
F1 = Smax1 / Smin1
F2 = Smax2 / Smin2
F3 = Smax3 / Smin3

zcrit = scipy.stats.f.ppf(0.975, k1, k2)

t11 =[S2x1, S2x2, k1, k2, F1]
t12 =[S2x1, S2x3, k1, k2, F2]
t13 =[S2x2, S2x3, k1, k2, F3]

if F1 > zcrit:
    t21 = [F1, zcrit, 0]
else:
    t21 = [F1, zcrit, 1]

if F2 > zcrit:
    t22 = [F2, zcrit, 0]
else:
    t22 = [F2, zcrit, 1]

if F3 > zcrit:
    t23 = [F3, zcrit, 0]
else:
    t23 = [F3, zcrit, 1]

f.write(str("Задание 6"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t11]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t12]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t13]))
f.write('\n')

f.write(str("Таблица 2"))
f.write('\n')
f.write(str("(1,2)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t21]))
f.write('\n')
f.write(str("(1,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t22]))
f.write('\n')
f.write(str("(2,3)"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t23]))
f.write('\n')

"Задание 7"
alpha = 0.05
t1, pval = scipy.stats.bartlett (x1, x2, x3)
if pval < alpha:
    t1 = [pval, alpha, 0]
else:
    t1 = [pval, alpha, 1]

f.write(str("Задание 7"))
f.write('\n')

f.write(str("Таблица 1"))
f.write('\n')
f.write(str([round(elem, 5) for elem in t1]))
f.write('\n')
