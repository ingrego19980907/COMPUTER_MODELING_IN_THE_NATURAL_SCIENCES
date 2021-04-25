import numpy as np
import matplotlib.pyplot as plt

# Всі генератори базуються на лінійному конгруентному: задаємо основні параметри
z = 2 ** 32
m = 69069
n = 1

# Задаємо початкове значення
x0 = 12


# Створимо функцію, що буде генерувати n випадкових рівномірно розподілених чисел від 0 до 1
def Rand01(N):  # n - кількість генерованих чисел
    R = np.zeros(N)
    x = x0
    for i in range(N):
        R[i] = round(x / z, 4)
        x = (m * x + n) % z
    return R


print(Rand01(10))

# Рівномірний розподіл від 10 до 15
# Спосіб 1. Генеруємо числа від 0 до 1, а потім виконуємо над ними потрібні дії
A = Rand01(10)
a = 10
b = 15
B = a + (b - a) * A  # Формула для перетворення закону розподілу, з лекції
print(B)


# Рівномірний розподіл від 10 до 15
# Спосіб 2. Створимо функцію

def Randab(N, a, b):
    R = np.zeros(N)
    x = x0
    for i in range(N):
        R[i] = round(x / z, 4)
        x = (m * x + n) % z
    Rab = a + (b - a) * R  # Формула для перетворення закону розподілу , з лекції
    # R - масив чисел від 0 до 1
    # Rab - масив чисел від а до b
    return Rab


print(Randab(10, 10, 15))

# Генеруємо 1000 чисел
Vipadkovi = Randab(1000, 10, 15)

# Перевіримо, чи не виходять числа за задані діапазони
print('Maximum:', max(Vipadkovi))
print('Minimum:', min(Vipadkovi))


# Створює функцію густини ймовірності чисел у масиві R, якщо весь діапазон значень чисел поділити на N ділянок
def NumR(a, b, R):
    R = np.array(R)
    return (((a <= R) & (R < b)).sum())


def Distribution(R, N):
    Rmin = min(R)
    Rmax = max(R)
    h = round((Rmax - Rmin) / N, 2)
    Levels = [round(Rmin + h * i + h / 2, 2) for i in range(N)]
    Probab = [round(NumR(Levels[i] - h / 2, Levels[i] + h / 2, R) / len(R) / h, 3) for i in range(N)]
    return (Levels, Probab)


# Проаналізуємо отриману послідовність на 20 проміжках
D = Distribution(Vipadkovi, 20)
plt.plot(D[0], D[1])
plt.title('Continuous uniform distribution from 10 to 15')
plt.ylabel('Probability density function')
plt.xlabel('Random number')
plt.xlim(10, 15)
plt.ylim(0, 2)
plt.grid(True)

D = Distribution(Randab(100000, 10, 15), 20)
plt.plot(D[0], D[1])
plt.title('Continuous uniform distribution from 10 to 15')
plt.ylabel('Probability density function')
plt.xlabel('Random number')
plt.xlim(10, 15)
plt.ylim(0, 0.3)
plt.grid(True)
