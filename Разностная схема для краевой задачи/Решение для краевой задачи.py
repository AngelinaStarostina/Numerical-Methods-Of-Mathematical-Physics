def mu_2(t, h, tau):
    return (2-h) * t + h / 2 + tau - 2 * h ** 2 / 5

def phi(x, t, tau):
    return x**2 - 2*t - tau


def solve(h, tau):
    stab = []
    N1 = int(1 / h)
    N2 = int(1 / tau)
    sigma = 1 / 2 - h ** 2 / (5 * tau)
    A_value = -tau * sigma / (h ** 2)
    C_value = 1 + tau * 2 * sigma / (h ** 2)
    B_value = -tau * sigma / (h ** 2)
    if abs(C_value) >= abs(B_value) + abs(A_value):
        stab.append(True)
    if abs(1 + h ** 2 / (2 * tau)) > 1:
        stab.append(True)
    for i in stab:
        if i == False:
            print('Достаточное условие устойчивости метода прогонки выполнено')
            return -1
    print('Достаточное условие устойчивости метода прогонки выполнено')

    y = []
    y.append([0 for j in range(N1 + 1)])
    for i in range(N2):
        A = [A_value for j in range(N1 + 1)]
        B = [B_value for j in range(N1 + 1)]
        C = [C_value for j in range(N1 + 1)]
        F = [0 for j in range(N1 + 1)]

        C[0] = 1
        A[0] = 0
        B[0] = 0
        for j in range(1, N1):
            F[j] = y[i][j] + ((1 - sigma) * tau / h ** 2) * (y[i][j - 1] - 2 * y[i][j] + y[i][j + 1]) + tau * phi(j * h, i * tau, tau)
        A[N1] = 1
        C[N1] = -1 - h**2 / (2 * tau * sigma)
        F[N1] = -(h * mu_2(tau * (i-1), h, tau) / sigma - (1-sigma)/sigma*(y[i-1][N1]-y[i-1][N1-1])+h**2/(2*tau*sigma)*y[i-1][N1])
        y.append(sweep_method(A, C, B, F))
    return y


def sweep_method(A, C, B, F):
    n = len(A)
    x = [0] * n
    for i in range(1, n):
        coef = A[i] / C[i - 1]
        C[i] = C[i] - coef * B[i - 1]
        F[i] = F[i] - coef * F[i - 1]
    x[n-1] = F[n - 1] / C[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = (F[i] - B[i] * x[i + 1]) / C[i]
    return x


y1 = solve(0.05, 0.0015)
N1 = int(1 / 0.05) + 1
N2 = int(1 /0.0015) + 1
for i in range(0, N2, 29):
    for j in range(0, N1, 2):
        print('{:.3e}'.format(y1[i][j]), end=" ")
    print()

print()
y2 = solve(0.05, 0.000015)
print("\nПогрешность:")
for i in range(0, N2, 29):
    for j in range(0, N1, 2):
        print('{:.3e}'.format(abs(y1[i][j] - y2[i*100][j])), end=" ")
    print()

