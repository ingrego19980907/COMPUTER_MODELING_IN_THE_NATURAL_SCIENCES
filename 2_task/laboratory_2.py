from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import math
import time


# Variant-18
# Двовимірна задача теплопровідності з джерелом тепла з граничними умовами
# Граничні умови задають сталу температуру на краях плівки

def c(j, i):
    if 0 <= j < 2:
        c = 0.05
    elif 2 <= j < 4:
        c = 0.02
    elif 4 <= j < 6:
        c = 0.01
    elif 6 <= j < 8:
        c = 0.07
    else:
        c = 0.1
    return c



# Задаємо кількість точок по часу
t0 = 0
tf = 20
h = 0.01
nt = math.floor((tf - t0) / h)
nt_ad = 400
nt1 = math.floor(nt / nt_ad)

# Задаємо кількість точок по координатах
L = 10
l = 0.25
nx = math.floor(L / l)
ny = nx
x = np.array([round(l * i, 3) for i in range(nx + 1)])
y = np.array([round(l * i, 3) for i in range(ny + 1)])

# Задаємо заготовку додаткового масиву, який будемо заповнювати
u1 = [[[0 for i in range(nx + 1)] for j in range(ny + 1)] for k in
      range(nt_ad + 1)]  # У такому масиві нумерація елементів йтиме як kji

# Задаємо початкові та крайові значення
for i in range(0, nx + 1):
    for j in range(0, ny + 1):
        u1[0][j][i] = 10
for k in range(nt_ad + 1):
    for i in [0, nx]:
        for j in range(0, ny + 1):
            u1[k][j][i] = 10
    for j in [0, ny]:
        for i in range(0, nx + 1):
            u1[k][j][i] = 10

# На даний момент в нульовому елементі додаткового масиву містяться початкові умови,які мають бути в основному масиві
u = [u1[0]]


# Задаємо джерело
def function(k, j, i):
    Q = 2 * h / c(j, i)
    return (Q)


t1 = time.process_time_ns()

# Задаємо великий цикл
for m in range(1):
    # Проводимо обчислення додаткового масиву
    for k in range(nt_ad):
        for j in range(1, ny):
            for i in range(1, nx):
                u1[k + 1][j][i] = round(h / c(j, i) / l / l * (u1[k][j][i + 1] - 2 * u1[k][j][i]
                                                               + u1[k][j][i - 1] + u1[k][j + 1][i] - 2 *
                                                                u1[k][j][i] + u1[k][j - 1][i])
                                                               + u1[k][j][i] + function(k, j, i), 5)
    # Записуємо 200-тий елемент в додатковому масиві до основного масиву
    u.append([[u1[-1][j][i] for i in range(nx + 1)] for j in range(ny + 1)])
    # Задаємо нове початкове значення для додаткового масиву
    # u1[0]=u1[-1]
    # І продовжуємо цикл
    print(m)
t2 = time.process_time_ns()

print('Час обчислень', round((t2 - t1) / 1000000000, 2), 'с')

N = 100
zcont = [10 + 50 * i / N for i in range(N + 1)]
zcol = [(1.0 * i / N, 0.9 - 0.9 * i / N, 1.0 - 1.0 * i / N) for i in range(N + 1)]

fig = plt.figure(figsize=(10, 10))

plt.contourf(x, y, u1[-1], zcont, colors=zcol)
plt.colorbar()
