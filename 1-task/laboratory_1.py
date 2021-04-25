from analysis_methods import Euler, EulerImp, Runge, Adams2E, Adams2I, Tochnist1, Tochnist1Ad2I, \
    Stiykist, StiykistAdams2I, Euler2, EulerCromer2, Midpoint2, Verlet2, VeloVerlet2, RelEr, RelErV, \
    Runge2, ProgonkaX, ProgonkaY, Runge_2
import math
import matplotlib.pyplot as plt

H = [0.01, 0.05, 0.1, 0.3, 0.5]
TF = 10
X0 = 1


def f1(t, x):
    return 0.4 * t - 3 * math.sqrt(t) + math.exp(0.3 * t)  # Функція, задана дифрівнянням варіант 18


def solution(t):
    return -4 / 3 + 10 / 3 * math.exp(0.3 * t) - 2 * t ** (3 / 2) + 0.2 * t * t


def f1ImpEu(t, x, h):
    return x + (0.4 * t - 3 * math.sqrt(t) + math.exp(0.3 * t)) * h


def f1ImpAd(t, x, x1, h):
    return (x1 + 5 / 3 * h * (t + 2 * h) + 2 / 3 * h * f1((t + h), x1) - h / 12 * f1(t, x)) / (
            1 + 1.25 * h - 2.5 / 12 * h * (t + 2 * h))


# Завдання 1.

# Аналіз стійкості явного методу Ейлера
Stiykist(Euler, f1, solution, H, 0, TF, X0)

# Аналіз точності явного методу Ейлера
Tochnist1(Euler, f1, solution, H, 0, TF, X0)

# Аналіз стійкості неявного методу Ейлера
Stiykist(EulerImp, f1ImpEu, solution, H, 0, TF, X0)

# Аналіз точності неявного методу Ейлера
Tochnist1(EulerImp, f1ImpEu, solution, H, 0, TF, X0)

# Аналіз стійкості явного методу Рунге
Stiykist(Runge, f1, solution, H, 0, TF, X0)

# Аналіз точності явного методу Рунге
Tochnist1(Runge, f1, solution, H, 0, TF, X0)

# Аналіз стійкості явного методу Адамса
Stiykist(Adams2E, f1, solution, H, 0, TF, X0)

# Аналіз точності явного методу Адамса
Tochnist1(Adams2E, f1, solution, H, 0, TF, X0)

# Аналіз стійкості неявного методу Адамса
StiykistAdams2I(Adams2I, f1ImpAd, f1, solution, H, 0, TF, X0)

# Аналіз точності неявного методу Адамса
Tochnist1Ad2I(Adams2I, f1ImpAd, f1, solution, H, 0, TF, X0)

# Завдання 2

t0 = 0
tf = 3
x0 = 0
x10 = 0


def fzavd2(t, x, x1):
    return 7 * x1 + 13 + 10 * math.cos(5 * t)
    # Функція, задана дифрівнянням у 2 завданні, яка визначає прискорення


def solution2zavd(t):
    return round((-3367 * t + 726 * math.exp(7 * t) - 343 * math.sin(5 * t) - 245 * math.cos(5 * t) - 481) / 1813)


def Verlet(t, xi, xi1, x1i, h):
    f = ((2 * xi1 - xi + h * h * (3.5 * xi / h + 13 + 10 * math.cos(5 * t))) / (1 - 2.5 * h))
    # Функція, виведена з дифрівняння 2 порядку для методу Верле
    return f


def VeloVerlet(t, xi, xi1, x1i, h):
    f = (x1i + 0.5 * h * (fzavd2(t, xi, x1i) + 13 + 10 * math.cos(5 * (t + h))) / (1 - 2.5 * h))
    # Функція, виведена з дифрівняння 2 порядку для методу Верле
    return f


def f1ImpAd(t, x, x1, h):
    return (x1 + 5 / 3 * h * (t + 2 * h) + 2 / 3 * h * f1((t + h), x1) - h * 12 * f1(t, x)) / (
            1 + 1.25 * h - 0.25 / 12 * h * (t + 2 * h))
    # Функція, виведена з дифрівняння 1 для неявного методу Адамса


(T, X, X1, Opis) = Euler2(fzavd2, 20, t0, tf, x0, x10)
S = [solution2zavd(T[i]) for i in range(len(T))]
print('X= ', X)
print('EulerCromer2: ', EulerCromer2(fzavd2, 20, t0, tf, x0, x10)[1])
print('Midpoint2: ', Midpoint2(fzavd2, 20, t0, tf, x0, x10)[1])
print('Verlet2: ', Verlet2(fzavd2, Verlet, 20, t0, tf, x0, x10)[1])
print('VeloVerlet2: ', VeloVerlet2(fzavd2, VeloVerlet, 20, t0, tf, x0, x10)[1])
print('S= ', S)
print('RelEr: ', RelEr(Midpoint2, fzavd2, solution2zavd, 20, t0, tf, x0, x10))
print('RelErV: ', RelErV(Verlet2, fzavd2, Verlet, solution2zavd, 20, t0, tf, x0, x10))

n0 = 20
for i in range(10):
    n = n0 * (2 ** i)
    E = RelEr(Midpoint2, fzavd2, solution2zavd, n, t0, tf, x0, x10)
    print('n=', n, 'Error = ', E)
    if E <= 0.1:
        break

n0 = 20
for i in range(8):
    n = n0 * (2 ** i)
    E = RelErV(Verlet2, fzavd2, Verlet, solution2zavd, n, t0, tf, x0, x10)
    print('n=', n, 'Error = ', E)
    if E <= 0.1:
        break


# Завдання 2: розв'язок системи
def fzavd2syst(t, x, x1):
    fun_1 = 7 * x1 + 13 + 10 * math.cos(5 * t)
    fun_2 = x1
    return fun_1, fun_2


n0 = 20
for i in range(10):
    n = n0 * (2 ** i)
    E = RelEr(Runge2, fzavd2syst, solution2zavd, n, t0, tf, x0, x10)
    print('n=', n, 'Error = ', E)
    if E <= 0.01:
        break

# Завдання 3

# Задаємо умови: варіант 18
k = 0.13
l = 52
xf = l
x0 = 0
y0 = 0
yf = 0
t0 = 0
tf = 6.7
n = 20

(T, X) = ProgonkaX(t0, tf, x0, xf, n, k)
Y = ProgonkaY(t0, tf, y0, yf, n, k)[1]

plt.plot(X, Y)


# Метод стрільби
# Розв'яжемо як систему методом Рунге-Кути

def fzavd3X(t, x, x1):
    f1 = (-k * x1)  # Функція, задана дифрівнянням у 2 завданні, яка визначає прискорення
    f2 = x1  # Функція для швидкості, отримана з представлення рівняння руху у вигляді системи
    return (f1, f2)


def fzavd3Y(t, y, y1):
    f1 = (-k * y1 - 9.8)  # Функція, задана дифрівнянням у 2 завданні, яка визначає прискорення
    f2 = y1  # Функція для швидкості, отримана з представлення рівняння руху у вигляді системи
    return (f1, f2)


# Задаємо умови: варіант 18
k = 0.13
l = 52
xf = l
x0 = 0
y0 = 0
yf = 0
t0 = 0
tf = 6.7
n = 20

# Підбираємо початкову швидкість - перший вистріл
x10 = 5
y10 = 10

XS = Runge2(fzavd3X, n, t0, tf, x0, x10)[1]
YS = Runge2(fzavd3Y, n, t0, tf, y0, y10)[1]

print(XS)
print(YS)

# Розв'яжемо для Х
# Метод звуження інтервалу вдвічі

x101 = 0.5  # При цьому значенні недолітаємо, швидкість пристрілки
x102 = 20  # При цьому значенні перелітаємо

# Задамо цикл на звуження проміжку вдвічі поки не прилетимо куди слід
d = 0.001  # Допустима похибка координати в кінці руху
for i in range(30):
    x103 = math.tan((math.atan(x101) + math.atan(x102)) / 2)
    C = Runge2(fzavd3X, n, t0, tf, x0, x103)[1][-1]
    # print('Швидкість',x103,'кінцева координата',C)
    if abs((C - xf)) < d:
        x10 = x103
        print(x103)
        break
    elif (C - xf) < 0:
        x101 = x103
    else:
        x102 = x103

# Розв'яжемо для Y
# Метод звуження інтервалу вдвічі

y101 = 10  # При цьому значенні недолітаємо, швидкість пристрілки
y102 = 20  # При цьому значенні перелітаємо

# Задамо цикл на звуження проміжку вдвічі поки не прилетимо куди слід
d = 0.001  # Допустима похибка координати в кінці руху
for i in range(30):
    y103 = math.tan((math.atan(y101) + math.atan(y102)) / 2)
    C = Runge2(fzavd3Y, n, t0, tf, y0, y103)[1][-1]
    # print('Швидкість',y103,'кінцева координата',C)
    if abs((C - yf)) < d:
        y10 = y103
        print(y103)
        break
    elif (C - yf) < 0:
        y101 = y103
    else:
        y102 = y103

# Вводимо знайдені початкові швидкості методом звуження інтервалу
x10 = 15.001481695995
y10 = 1.320018948282612

XS = Runge2(fzavd3X, n, t0, tf, x0, x10)[1]
YS = Runge2(fzavd3Y, n, t0, tf, y0, y10)[1]

print(XS)
print(YS)

plt.plot(XS, YS)

# Розв'яжемо для Х
# Метод для лінійних задач

# Робимо дві пристрілки: будь-які
x101 = 200
x102 = 300

XS1 = Runge2(fzavd3X, n, t0, tf, x0, x101)[1]
XS2 = Runge2(fzavd3X, n, t0, tf, x0, x102)[1]

B1 = XS1[-1]
B2 = XS2[-1]  # Відповідно до формул лекція 3-4 слайд 45

# Створюємо таблицю розв'язку як суму двох пристрілочних розв'язків
XStr = [((xf - B2) * XS1[i] + (B1 - xf) * XS2[i]) / (B1 - B2) for i in range(len(XS1))]

# Розв'яжемо для Y
# Метод для лінійних задач

# Робимо дві пристрілки
y101 = 10
y102 = 20

YS1 = Runge2(fzavd3Y, n, t0, tf, y0, y101)[1]
YS2 = Runge2(fzavd3Y, n, t0, tf, y0, y102)[1]

B1 = YS1[-1]
B2 = YS2[-1]  # Відповідно до формул лекція 3-4 слайд 45

# Створюємо таблицю розв'язку як суму двох пристрілочних розв'язків
YStr = [((yf - B2) * YS1[i] + (B1 - yf) * YS2[i]) / (B1 - B2) for i in range(len(YS1))]

plt.plot(XStr, YStr)


# Завдання 4

# Задаємо систему рівнянь

def func4zavd(t, x, y):
    funcX = a * (1 - x / k1) * x - b / (k2 - x) * x * y
    funcY = (c - d * y / x) * y
    return (funcX, funcY)


# Задаємо умови задачі варіант 31
x0 = 42
y0 = 3
a = 1.6
b = 3
c = 0.08
d = 0.05
k1 = 9.1
k2 = 1.9

# Будуємо графіки розв'язків
(T, X, Y, opys) = Runge_2(func4zavd, 2000, 0, 300, x0, y0)

print(X)
print(Y)
plt.plot(T, X)
plt.plot(T, Y)
