import matplotlib.pyplot as plt
import numpy as np
import math


# Явний метод Ейлера
def Euler(function, h, t0, tf, x0):
    n = math.ceil((tf - t0) / h)  # обчислюємо кількість кроків
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):
        x = x + [round(x[i] + h * function(t[i], x[i]), 3)]
    return t, x, (h, t0, tf, x0, 'Explicit Euler')  # Результат, який отримаємо при виклику функції


# Неявний метод Ейлера
def EulerImp(function_implicit, h, t0, tf, x0):
    n = math.ceil((tf - t0) / h)  # обчислюємо кількість кроків
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):
        x = x + [round(function_implicit(t[i + 1], x[i], h), 3)]
    return t, x, (h, t0, tf, x0, 'Implicit Euler')


# Метод Рунге-Кути
def Runge(function, h, t0, tf, x0):
    n = math.ceil((tf - t0) / h)
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    # x1=[] Рядок формує заготовку для таблиці швидкостей зростання функції - 1 похідна
    for i in range(n):
        k1 = function(t[i], x[i])  # Обчислюємо к1 за м. Рунге-Кути 4 порядку
        k2 = function(t[i] + h / 2, x[i] + h * k1 / 2)
        k3 = function(t[i] + h / 2, x[i] + h * k2 / 2)
        k4 = function(t[i] + h, x[i] + h * k3)
        x = x + [round(x[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 3)]
    return t, x, (h, t0, tf, x0, 'Runge 4 order')


# Метод Адамса
def Adams2E(function, h, t0, tf, x0):
    n = math.ceil((tf - t0) / h)
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    # x1=[] Рядок формує заготовку для таблиці швидкостей зростання функції - 1 похідна
    # Перший крок робимо за Рунге
    i = 0
    k1 = function(t[i], x[i])  # Обчислюємо к1 за м. Рунге-Кути 4 порядку
    k2 = function(t[i] + h / 2, x[i] + h * k1 / 2)
    k3 = function(t[i] + h / 2, x[i] + h * k2 / 2)
    k4 = function(t[i] + h, x[i] + h * k3)
    x = x + [round(x[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 3)]
    for i in range(n - 1):  # Заходимо в цикл, що складається з (n-1) кроків
        x = x + [round(x[i + 1] + 3 / 2 * h * function(t[i + 1], x[i + 1]) - h / 2 * function(t[i], x[i]), 3)]
    return t, x, (h, t0, tf, x0, 'Adams 2-points Explicit')


# Неявний метод Адамса
def Adams2I(functionImplicitAdams, function, h, t0, tf, x0):
    n = math.ceil((tf - t0) / h)
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    # Перший крок робимо за Рунге
    i = 0
    k1 = function(t[i], x[i])  # Обчислюємо к1 за м. Рунге-Кути 4 порядку
    k2 = function(t[i] + h / 2, x[i] + h * k1 / 2)
    k3 = function(t[i] + h / 2, x[i] + h * k2 / 2)
    k4 = function(t[i] + h, x[i] + h * k3)
    x = x + [round(x[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), 3)]
    for i in range(n - 1):
        x = x + [round(functionImplicitAdams(t[i], x[i], x[i + 1], h), 3)]
    return t, x, (h, t0, tf, x0, 'Adams 2-points Implicit')


def Tochnist1(Method, function, solution, H, t0, tf, x0):
    # задаємо назву методу, буде розв'язувати задане рівняння при даній початковій умові з різними кроками з масиву Н,
    # додатково треба буде задати аналітичний розв'язок
    Error = []
    for i in range(len(H)):
        # Входимо в цикл, кількість повторень - кількість елементів в масиві кроків, які хочемо проаналізувати
        (T, X, Opis) = Method(function, H[i], t0, tf, x0)  # З і-м кроком розв'язуємо рівняння заданим методом
        Sol = [round(solution(T[i]), 3) for i in range(len(T))]  # Таблиця значень функції аналітичного розв'язку
        Error = Error + [round(math.sqrt(np.sum([(Sol[i] - X[i]) ** 2 for i in range(len(T))]) / (len(T))), 5)]
        # Додаємо до таблиці похибок новобчислений елемент, що міститься в квадратних дужках
        # Округлений до 5 знаків після коми орінь з відношення суми елементів масиву
        # ([(Sol[i]-X[i])**2 for i in range(len(T))]), що складається з квадратів різниць точних та чисельних значень,
        # до кількості елементів
    fig, a = plt.subplots()
    a.plot(H, Error)
    plt.title('Dependence of the error of ' + Opis[4] + ' method on time step')
    plt.ylabel('Error')
    plt.xlabel('Time step')
    plt.grid(True)
    return H, Error, (t0, tf, x0, Opis[4])


def Tochnist1Ad2I(Method, functionImplicitAdams, function, solution, H, t0, tf, x0):
    # задаємо назву методу, буде розв'язувати задане рівняння при даній початковій умові з різними кроками з масиву Н,
    # додатково треба буде задати аналітичний розв'язок
    Error = []  # Створюємо порожній масив, як заготовку для значень похибки
    for i in range(len(H)):
        (T, X, Opis) = Method(functionImplicitAdams, function, H[i], t0, tf, x0)
        Sol = [round(solution(T[i]), 3) for i in range(len(T))]
        Error = Error + [round(math.sqrt(np.sum([(Sol[i] - X[i]) ** 2 for i in range(len(T))]) / (len(T))), 5)]
        # Додаємо до таблиці похибок новобчислений елемент, що міститься в квадратних дужках
        # Округлений до 5 знаків після коми орінь з відношення суми елементів масиву
        # ([(Sol[i]-X[i])**2 for i in range(len(T))]), що складається з квадратів різниць точних та чисельних значень,
        # до кількості елементів
    fig, a = plt.subplots()
    a.plot(H, Error)
    plt.title('Dependence of the error of ' + Opis[4] + ' method on time step')
    plt.ylabel('Error')
    plt.xlabel('Time step')
    plt.grid(True)
    return H, Error, (t0, tf, x0, Opis[4])


def Stiykist(Method, function, solution, H, t0, tf, x0):
    plt.ylabel('x')  # Підпис горизонтальної осі
    plt.xlabel('t')  # Підпис вертикальної осі
    plt.grid(True)  # Увімкнення сітки
    plt.xlim(t0, tf)  # Задає межі по горизонталі
    for i in range(len(H)):
        (T, X, Opis) = Method(function, H[i], t0, tf, x0)  # З і-м кроком розв'язуємо рівняння заданим методом
        if i == 0:
            Sol = [round(solution(T[i]), 3) for i in range(len(T))]  # Таблиця значень функції аналітичного розв'язку
            plt.plot(T, Sol, label='Analytical')  # Будуємо графік аналітичного розв'язку
        plt.plot(T, X, label='h=' + str(H[i]))  # Будуємо графік розв'язку
    plt.title(Opis[4] + ' with varying step')  # Підписуємо графік
    plt.legend()  # Вмикаємо легенду
    return ()


def StiykistAdams2I(Method, functionImplicitAdams, function, solution, H, t0, tf, x0):
    plt.ylabel('x')  # Підпис горизонтальної осі
    plt.xlabel('t')  # Підпис вертикальної осі
    plt.grid(True)  # Увімкнення сітки
    # plt.ylim(0,1.0) Рядок задавав би межі по вертикалі, ми зробили їх автоматичними
    plt.xlim(t0, tf)  # Задає межі по горизонталі
    for i in range(len(H)):
        (T, X, Opis) = Method(functionImplicitAdams, function, H[i], t0, tf, x0)
        if i == 0:
            Sol = [round(solution(T[i]), 3) for i in range(len(T))]  # Таблиця значень функції аналітичного розв'язку
            plt.plot(T, Sol, label='Analytical')  # Будуємо графік аналітичного розв'язку
        plt.plot(T, X, label='h=' + str(H[i]))  # Будуємо графік розв'язку
    plt.title(Opis[4] + ' with varying step')  # Підписуємо графік
    plt.legend()  # Вмикаємо легенду
    return ()


# Явний метод Ейлера для 2 порядку
def Euler2(function, n, t0, tf, x0, x10):
    # function - функція, задана дифрівнянням, яка визначає швидкість приросту координати
    h = round((tf - t0) / n, 3)  # обчислюємо крок
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):  # Заходимо в цикл, що складається з n кроків
        x1 = x1 + [round(x1[i] + function(t[i], x[i], x1[i]) * h, 3)]
        x = x + [round(x[i] + h * x1[i], 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Explicit Euler')


# Напівнеявний метод Ейлера для 2 порядку (метод Ейлера-Кромера)
def EulerCromer2(function, n, t0, tf, x0, x10):
    h = round((tf - t0) / n, 3)  # обчислюємо крок
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):  # Заходимо в цикл, що складається з n кроків
        x1 = x1 + [round(x1[i] + function(t[i], x[i], x1[i]) * h, 3)]
        x = x + [round(x[i] + h * x1[i + 1], 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Euler Cromer')


# Метод серединної точки для 2 порядку
def Midpoint2(function, n, t0, tf, x0, x10):
    h = round((tf - t0) / n, 3)  # обчислюємо крок
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):
        x1 += [round(x1[i] + function(t[i], x[i], x1[i]) * h, 3)]
        x = x + [round(x[i] + h * (x1[i] + x1[i + 1]) / 2, 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Midpoint')


# Метод Верле
def Verlet2(function, functionVerlet, n, t0, tf, x0, x10):
    h = round((tf - t0) / n, 3)  # обчислюємо крок
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    i = 0  # Перший крок - методом серединної точки
    x1 = x1 + [round(x1[i] + function(t[i], x[i], x1[i]) * h, 3)]
    x = x + [round(x[i] + h * (x1[i] + x1[i + 1]) / 2, 3)]
    for i in range(n - 1):  # Заходимо в цикл, що складається з (n-1) кроків
        x = x + [round(functionVerlet(t[i], x[i], x[i + 1], x1[i], h), 3)]
        if i != 0:
            x1 = x1 + [round((x[i + 2] - x[i]) / 2 / h, 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Verlet Integration')


# Метод Верле швидкісний
def VeloVerlet2(function, funcVeloVerlet, n, t0, tf, x0, x10):
    h = round((tf - t0) / n, 3)  # обчислюємо крок
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):  # Заходимо в цикл, що складається з (n) кроків
        x = x + [round(x[i] + x1[i] * h + function(t[i], x[i], x1[i]) * h * h / 2, 3)]
        x1 = x1 + [round(funcVeloVerlet(t[i], x[i], x[i + 1], x1[i], h), 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Verlet Integration')


# Функція оцінки похибки
def RelEr(Method, function, solution, n, t0, tf, x0,
          x10):  # задаємо назву методу, буде розв'язувати задане рівняння при даній початковій умові з різними кроками з масиву Н, додатково треба буде задати аналітичний розв'язок
    (T, X, X1, Opis) = Method(function, n, t0, tf, x0, x10)  # З і-м кроком розв'язуємо рівняння заданим методом
    sol = [round(solution(T[i]), 3) for i in range(len(T))]  # Таблиця значень функції аналітичного розв'язку
    relError = round(
        math.sqrt(np.sum([((sol[i] - X[i]) / sol[i] if sol[i] != 0 else 0) ** 2 for i in range(len(T))]) / (len(T))), 5)
    return relError


# Функція оцінки похибки для методу Верле
def RelErV(Method, function, functionVerlet, solution, n, t0, tf, x0, x10):
    (T, X, X1, Opis) = Method(function, functionVerlet, n, t0, tf, x0, x10)
    sol = [round(solution(T[i]), 3) for i in range(len(T))]
    relError = round(math.sqrt(np.sum([((sol[i] - X[i]) / sol[i] if sol[i] != 0 else 0) ** 2 for i in range(len(T))]
                                      ) / (len(T))), 5)
    return relError


# Метод Рунге-Кути для системи з 2 рівнянь
def Runge2(function, n, t0, tf, x0, x10):
    h = round((tf - t0) / n, 3)
    t = [round(t0 + i * h, 3) for i in range(n + 1)]  # Створюємо таблицю дискретного часу
    x = [x0]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    x1 = [x10]  # Створюємо масив, який поки що містить 1 елемент (нульовий елемент - початкова умова)
    for i in range(n):  # Заходимо в цикл, що складається з n кроків
        k1x1 = function(t[i], x[i], x1[i])[0]  # Обчислюємо к1 для швидкості за м. Рунге-Кути 4 порядку
        k1x = function(t[i], x[i], x1[i])[1]  # Обчислюємо к1 для координати за м. Рунге-Кути 4 порядку
        (k2x1, k2x) = function(t[i] + h / 2, x[i] + h * k1x / 2, x1[i] + h * k1x1 / 2)
        (k3x1, k3x) = function(t[i] + h / 2, x[i] + h * k2x / 2, x1[i] + h * k2x1 / 2)
        (k4x1, k4x) = function(t[i] + h, x[i] + h * k3x, x1[i] + h * k3x1)
        x = x + [round(x[i] + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x), 3)]
        x1 = x1 + [round(x1[i] + h / 6 * (k1x1 + 2 * k2x1 + 2 * k3x1 + k4x1), 3)]
    return t, x, x1, (h, t0, tf, x0, x10, 'Runge 4 order')
