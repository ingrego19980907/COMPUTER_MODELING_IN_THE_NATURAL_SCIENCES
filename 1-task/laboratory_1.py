from analysis_methods import Euler, EulerImp, Runge, Adams2E, Adams2I, Tochnist1, Tochnist1Ad2I, \
    Stiykist, StiykistAdams2I, Euler2, EulerCromer2, Midpoint2, Verlet2, VeloVerlet2, RelEr, RelErV, Runge2
import math

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
