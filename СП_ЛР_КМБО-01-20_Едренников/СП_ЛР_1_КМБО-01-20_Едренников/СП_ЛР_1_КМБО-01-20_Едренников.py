import numpy
import random
from decimal import Decimal

f = open('answer.txt', 'r+')

r = [[-2, -1, -1, -1, -1, -1, -1, -1, -1, -1],
     [0, -1, 0, 0.332, 0, 0, 0, 0, 0, 0],
     [0.332, 0.332, -1, 0, 0, 0, 0, 0, 0, 0],
     [0, 0, 0, -1, 0.332, 0.332, 0, 0, 0, 0],
     [0, 0, 0.332, 0, -1, 0.069, 0, 0, 0, 0],
     [0.069, 0.069, 0.069, 0.069, 0.069, -1, 0, 0, 0, 0],
     [0, 0, 0, 0, 0, 0, -1, 0.332, 0.332, 0.332],
     [0.599, 0, 0, 0, 0, 0, 0.332, -1, 0.069, 0.069],
     [0, 0.599, 0.599, 0, 0, 0, 0.069, 0.069, -1, 0.599],
     [0, 0, 0.599, 0.599, 0.599, 0.599, 0.599, 0.599, 0.599, -1],
     ]

ro = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]

rA = numpy.linalg.solve(r, ro)

print(rA)
print(sum(rA))

P = [[0, 0, 0.332, 0, 0.069, 0, 0, 0.599, 0, 0],
     [0, 0, 0.332, 0, 0, 0.069, 0, 0, 0.599, 0],
     [0, 0.332, 0, 0, 0, 0.069, 0, 0, 0.599, 0],
     [0, 0, 0, 0, 0.332, 0.069, 0, 0, 0, 0.599],
     [0, 0, 0, 0.332, 0, 0.069, 0, 0, 0, 0.599],
     [0, 0, 0, 0.332, 0.069, 0, 0, 0, 0, 0.599],
     [0, 0, 0, 0, 0, 0, 0, 0.332, 0.069, 0.599],
     [0, 0, 0, 0, 0, 0, 0.332, 0, 0.069, 0.599],
     [0, 0, 0, 0, 0, 0, 0.332, 0.069, 0, 0.599],
     [0, 0, 0, 0, 0, 0, 0.332, 0.069, 0.599, 0]]

pk = P
delt = 1
K = 1
deltak = [0]*10
while delt > 0.001:
    delt = 0
    pk = numpy.dot(pk, pk)
    for i in range(0, 10):
        for j in range(0, 10):
            delt = max((abs(pk[i][j] - rA[j])), delt)
    K += 1

for i in range(0, 10):
        for j in range(0, 10):
            deltak[i] = Decimal(max((abs(pk[i][j] - rA[j])), deltak[i]))

v = numpy.zeros((10, 10), dtype=float, order='C')
vn = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
for i in range(0, 10):
    t = i
    for j in range(0, 100):
        k = random.choices(vn, P[t])
        v[i][k[0] - 1] += 1
        t = k[0]-1

v = v / 100
vr = numpy.zeros((10, 10))

for i in range(0, 10):
    for j in range(0, 10):
        vr[i][j] = abs(v[i][j] - rA[j])

f.write(str(numpy.around(rA, 6)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(K))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(numpy.around(pk, 6)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(deltak))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(numpy.around(v, 6)))
f.write('\n')
f.write('\n')
f.write('\n')
f.write(str(numpy.around(vr, 6)))
f.write('\n')
