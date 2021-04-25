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
