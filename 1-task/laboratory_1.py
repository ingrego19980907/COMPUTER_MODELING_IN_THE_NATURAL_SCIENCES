from analysis_methods import Euler, EulerImp, Runge, Adams2E, Adams2I, Tochnist1, Tochnist1Ad2I, \
    Stiykist, StiykistAdams2I
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
Stiykist(EulerImp, f1, solution, H, 0, TF, X0)

# Аналіз точності неявного методу Ейлера
Tochnist1(EulerImp, f1, solution, H, 0, TF, X0)

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

