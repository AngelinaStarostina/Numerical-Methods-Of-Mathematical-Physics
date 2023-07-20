import math
import numpy as np

h = 0.1
N = int(1/h)
k0 = 2
g0 = 2
k1 = 0
g1 = -3 * math.sin(1)
y = np.zeros(N+1)
alpha = np.zeros(N)
beta = np.zeros(N)
a = np.zeros(N+1)

def x(num):
    return num * h

def f(i):
    if i == N:
        return 4 * math.cos(1 - h / 4) - (2 - h / 2) * math.sin(1 - h / 4)
    elif i == 0:
        return 4 * math.cos(h / 4) - (h / 2) * math.sin(h / 4)
    else:
        return 4 * math.cos(x(i)) - 2 * x(i) * math.sin(x(i))

def k(x):
    return 4 - x ** 2

def q(x):
    return x**2

def calculate_a():
    for i in range(N+1):
        a[i] = 4 - ((x(i) + x(i-1)) / 2) ** 2

def d(i):
    if i == 0:
        return h**2 / 16
    if i == N:
        return (1 - h / 4)**2
    else:
        return x(i)**2

def B(i):
    return a[i+1] / h**2


def A(i):
    return a[i] / h**2


def C(i):
    return (a[i+1] + a[i]) / h**2 + d(i)


def find_alpha():
    alpha[0] = a[0] / (h * (k0 + h/2 * d(0) + a[0]/h))
    for i in range(1, N):
        alpha[i] = B(i) / (C(i) - A(i) * alpha[i - 1])


def find_beta():
    beta[0] = (g0 + h/2*f(0)) / (k0 + h/2 * d(0) + a[0]/h)
    for i in range(1, N):
        beta[i] = (f(i) + beta[i-1] * A(i)) / (C(i) - A(i) * alpha[i-1])


def sweep_method():
    find_alpha()
    find_beta()
    k2_sweep = a[N] / (h * (k1 + h / 2 * d(N) + a[N] / h))
    nu2 = (g1 + h / 2 * f(N)) / (k1 + h / 2 * d(N) + a[N] / h)
    y[N] = (nu2 + k2_sweep * beta[N-1]) / (1 - alpha[N-1] * k2_sweep)
    for i in reversed(range(0, N)):
        y[i] = alpha[i] * y[i+1] + beta[i]


def test_stab():
    stab = []
    B_ = np.zeros(N)
    for i in range(1, N):
        B_[i] = B(i - 1)
    A_ = np.zeros(N)
    for i in range(1, N):
        A_[i] = A(i - 1)
    C_ = np.zeros(N)
    for i in range(1, N):
        C_[i] = C(i - 1)
    for i in range(N):
        stab.append(abs(C_[i]) >= abs(A_[i]) + abs(B_[i]))
    stab.append(alpha[0] <= 1)
    stab.append(a[N] / (h * (k1 + h / 2 * d(N) + a[N] / h)) <= 1)
    stab.append(abs(k(1) / (h * (k(1) / h + k1)) + abs(alpha[0]) < 2))
    for i in range(len(stab)):
        if stab[i] == False:
            print("Достаточное условие устойчивости метода прогонки не выполнено")
            return -1
    print("Достаточное условие устойчивости метода прогонки выполнено")
    return 0


test_stab()
calculate_a()
sweep_method()
print("Решение: ", y)

res = np.copy(y)
h = 0.001
N = int(1/h)
y = np.zeros(N+1)
alpha = np.zeros(N)
beta = np.zeros(N)
a = np.zeros(N+1)
calculate_a()
sweep_method()
res_exact = np.zeros_like(res)
for i in range(int(N/100)+1):
    res_exact[i] = y[i*100]
print("Решение при h = 0.001:", res_exact)

e = res_exact - res
print("Погрешность: ", e)