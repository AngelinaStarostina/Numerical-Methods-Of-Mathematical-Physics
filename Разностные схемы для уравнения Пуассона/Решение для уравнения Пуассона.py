import numpy as np

h1 = 0.1
h2 = 0.05
N1 = int(1 / h1)
N2 = int(1 / h2)
res = np.zeros((N2+1, N1+1))
eps = max(h1**2, h2**2)
dif = np.zeros((N2+1, N1+1))


def x(i):
    return i * h1


def y(j):
    return j * h2


def f(x, y):
    return 12 * (x**2 + y**2)


def fill_border():
    for i in range(N1+1):
        res[0, i] = x(i)**4
        res[N2, i] = 1 + x(i)**4

    for i in range(1, N2):
        res[i, 0] = y(i)**4
        res[i, N1] = 1 + y(i)**4


def count_max(max_matr):
    max = -1000
    for i in range(1, N2):
        for j in range(1, N1):
            if max_matr[i, j] > max:
                max = max_matr[i, j]
    return max


def iteration():
    y1 = np.copy(res)
    h1_2 = h1**2
    h2_2 = h2**2
    max_matr = np.zeros_like(y1)+11
    while count_max(max_matr) > eps:
        y2 = np.copy(y1)
        max_matr = np.zeros_like(y1)
        for j in reversed(range(1, N2)):
            for i in range(1, N1):
                y1[j][i] = 1 / (2/h1_2 + 2/h2_2) * ((y2[j+1, i] + y2[j-1, i])/h1_2 + (y2[j, i+1] + y2[j, i-1])/h2_2 + f(x(j), y(i)))
                max_matr[j, i] = abs(y1[j, i] - y2[j, i])
    return y1


fill_border()
res = iteration()
print('\n'.join(' '.join(map(str, ["{0:.11f}".format(x) for x in res[N2-i]])) for i in range(N2+1)))
print()


h1 /= 10
h2 /= 10
N1 = int(1 / h1)
N2 = int(1 / h2)
eps = max(h1**2, h2**2)
res2 = np.copy(res)
res = np.zeros((N2+1, N1+1))
fill_border()
res = iteration()
res3 = np.zeros_like(res2)
for i in range(int(N2/10)+1):
    for j in range(int(N1/10)+1):
        k = int(i*10)
        l = int(j*10)
        res3[i, j] = res[k, l]

print('\n'.join(' '.join(map(str, ["{0:.11f}".format(x) for x in res3[int(N2/10)-i]])) for i in range(int(N2/10)+1)))
print()

for i in range(int(N2/10)):
    for j in range(int(N1/10)):
        k = int(i*10)
        l = int(j*10)
        dif[i, j] = abs(res2[i, j]-res[k, l])
print('\n'.join(' '.join(map(str, ["{0:.11f}".format(x) for x in dif[int(N2/10-i)]])) for i in range(int(N2/10+1))))