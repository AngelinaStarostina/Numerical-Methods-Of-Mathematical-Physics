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
res_exact = [1.00000016, 0.99500433, 0.98006675, 0.95533666, 0.92106117, 0.87758274, 0.82533579, 0.76484237, 0.69670689, 0.62161015, 0.54030248]


def x(num):
    return num * h


def f(i):
    return 4 * math.cos(x(i)) - 2 * x(i) * math.sin(x(i))


def k(x):
    return 4 - x ** 2


def q(x):
    return x**2


def dk(x):
    return -2*x


def calculate_k0():
    return k0 * (1 - (h / 2) * (dk(0) / k(0))) + h / 2 * q(0)


def calculate_g0():
    return g0 * (1 - (h / 2) * (dk(0) / k(0))) + h / 2 * f(0)


def calculate_k1():
    return k1 * (1 + (h / 2) * (dk(1) / k(1))) + h / 2 * q(1)


def calculate_g1():
    return g1 * (1 + (h / 2) * (dk(1) / k(1))) + h / 2 * f(N)


def B(num):
    xi = x(num)
    return k(xi) / (h ** 2) + dk(xi) / (2 * h)


def A(num):
    xi = x(num)
    return k(xi) / (h ** 2) - dk(xi) / (2 * h)


def C(num):
    xi = x(num)
    return 2 * k(xi) / h**2 + q(xi)


def find_alpha():
    alpha[0] = k(0) / (k0 * h + k(0))
    for i in range(1, N):
        alpha[i] = B(i) / (C(i) - A(i) * alpha[i - 1])


def find_beta():
    beta[0] = g0 * h / (k0 * h + k(0))
    for i in range(1, N):
        beta[i] = (f(i) + beta[i-1] * A(i)) / (C(i) - A(i) * alpha[i-1])


def sweep_method():
    find_alpha()
    find_beta()
    k2_sweep = k(1) / (k1 * h + k(1))
    nu2 = g1 * h / (k1 * h + k(1))
    y[N] = (nu2 + k2_sweep * beta[N-1]) / (1 - alpha[N-1] * k2_sweep)
    for i in reversed(range(0, N)):
        y[i] = alpha[i] * y[i+1] + beta[i]


def print_1000():
    for i in range(int(N/100)):
        print(y[i*100])


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
    stab.append(k(1) / (h * (k(1) / h + k1)) <= 1)
    stab.append(abs(k(1) / (h * (k(1) / h + k1)) + abs(alpha[0]) < 2))
    for i in range(len(stab)):
        if stab[i] == False:
            print("Достаточное условие устойчивости метода прогонки не выполнено")
            return -1
    print("Достаточное условие устойчивости метода прогонки выполнено")
    return 0


test_stab()
g1 = calculate_g1()
g0 = calculate_g0()
k1 = calculate_k1()
k0 = calculate_k0()
sweep_method()
print("Решение: ", y)

e = res_exact - y
print("Погрешность: ", e)
