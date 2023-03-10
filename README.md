# Методы оптимизации. Лабораторная работа №1

**Создатели**: Бактурин Савелий, Сотников Максим, Вереня Андрей

## Постановка задачи и примерное решение

### Реализуйте градиентный спуск с постоянным шагом (learning rate)

#### Литература

http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BD%D0%BE%D0%B3%D0%BE_%D1%81%D0%BF%D1%83%D1%81%D0%BA%D0%B0

http://www.machinelearning.ru/wiki/images/b/b5/MO17_seminar5.pdf

#### Конспект

Что же такое градиентный спуск? Метод выглядит так. Возьмём $x_0$ - некое начальное приближение. Как его выбрать - тоже задача непростая, в зависимости от него мы получим разные результаты. Если мы примерно догадываемся, где минимум, можно попытаться выбрать точку близко к нему, например (хотя это не гарантирует нам быструю сходимость). Так вот, у нас есть начальное условие, и мы двигаем наши параметры по формуле

$$
  x_{k+1}=x_k+\alpha p_k
$$

При этом хороший вопрос~ - когда остановится. Есть разные соображения о том, когда останавливаться. Можно смотреть на условия минимума и как-то их проверять. Можно смотреть на Гессиан, если нам нужна высокая гарантия, что мы нашли минимум, но это делают достаточно редко, потому что это сложно. Но в любом случае это не делают единственным критерием. Вторым критерием делают изменение $x$. Если он меняется очень мало, с большой вероятностью мы в тупике.

Ещё нам часто хочется, чтобы функция убывала с ростом $k$ (это, опять же, не что-то универсальное, в стохастическом градиентном спуске это не выполняется, но тут хочется). Для этого надо наложить условия на $p_k$:
$$
  p_k^T\nabla f(x_k)<0
$$

Самое очевидное - в качестве $p_k$ взять $-\nabla f(x_k)$. Это, собственно, и есть градиентный спуск. Это не единственный вариант, что можно делать, но один из. Альтернативные пути: можно смотреть на предыдущие шаги каким-то образом (метод сопряжённых градиентов, например) или можно пользоваться Ньютоновскими шагами.

Вопрос: что такое $\alpha$ (в машинном обучении называют коэффициент обучения)? Простейший вариант~--- какая-то константа. Но также можно менять $\alpha$ в зависимости от номера шага (называется learning rate scheduling). Впрочем, такое не очень имеет смысл использовать в обычном градиентном спуске, помогает только в стохастическом.

Что ещё можно делать? Можно искать минимум вдоль направления. Тут в первую очередь возникают два вопроса: как и насколько это нам надо (типа, насколько это нужно, если это просто один из шагов алгоритма). На второй мы не ответим, а про первый немного поговорим. Для начала мы ищем интервал, где находится минимум. Как? Посмотрим на условие $p_k^T\nabla f(x_k)<0$. Это даёт направление, в котором функция убывает. Мы смотрим на направление, куда функция убывает, и идём туда с каким-то шагом до тех пор, пока функция убывает.  Скорее всего наша функция не убывает бесконечно, а значит когда начнёт возрастать, мы найдём три точки, где центральная меньше боковых, то есть наши отрезок, где лежит минимум (bracketing называется). Дальше можно использовать деление пополам (или какие-то его вариации типа золотого сечения или ещё чего). Можно взять точки, провести многочлен через них и взять минимум (это интерполяция). Можно совместить интерполяцию с дихотомией. Или можно использовать что-то, что использует производные.

### Реализуйте метод одномерного поиска (метод дихотомии, метод Фибоначчи, метод золотого сечения) и градиентный спуск на его основе.

#### Литература

http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4%D1%8B_%D0%B4%D0%B8%D1%85%D0%BE%D1%82%D0%BE%D0%BC%D0%B8%D0%B8

#### Что говорит Chat GPT, первый

Метод дихотомии:

1. Задать начальные значения a и b, где f(a) < 0 и f(b) > 0
2. Вычислить среднюю точку c = (a + b) / 2
3. Если f(c) близко к нулю, то c является корнем уравнения, и мы заканчиваем поиск
4. Если f(c) < 0, то корень находится в интервале (c, b), и мы переходим к шагу 2, задавая новые значения a и b
5. Если f(c) > 0, то корень находится в интервале (a, c), и мы переходим к шагу 2, задавая новые значения a и b

Градиентный спуск на основе метода дихотомии:

1. Задать начальное значение x0
2. Повторять до сходимости:
   - Вычислить градиент функции f(x) в точке xk
   - Найти направление спуска dk = -grad(f(xk))
   - Вычислить длину шага alpha методом дихотомии на отрезке [0, 1]
   - Обновить значение xk+1 = xk + alpha * dk
3. Вернуть xk+1 как приближенное решение

Примечание: для обеспечения сходимости градиентного спуска необходимо выбрать правильное начальное значение x0 и правильно настроить параметры метода (например, длину шага alpha).

#### Что говорит Chat GPT, второй

Sure, I can help you with that. 

For the one-dimensional search method (dichotomy method), the basic algorithm is as follows:
Step 1: Initialize a and b such that f(a) and f(b) have opposite signs.
Step 2: Let c = (a+b)/2, calculate f(c).
Step 3: If f(c) is zero, then stop, c is the root.
Step 4: If f(c) and f(a) have opposite signs, then set b = c, otherwise set a = c.
Step 5: Repeat steps 2-4 until desired accuracy is achieved. 

This algorithm essentially involves repeatedly bisecting an interval and determining which subinterval contains a root of the function. 

As for gradient descent, it is an iterative optimization algorithm for finding the minimum of a function. The idea is to take steps proportional to the negative of the gradient at the current point. The algorithm can be written as follows:

Step 1: Initialize x0 and choose a learning rate alpha.
Step 2: Repeat until convergence:
         - Compute the gradient of the function at x_k.
         - Update x_{k+1} = x_k - alpha * gradient at x_k.


The key idea behind gradient descent is that by iteratively taking steps in the direction of steepest descent, we can converge to the minimum of the function. 

Hope this helps! Let me know if you have any other questions.

### Проанализируйте траекторию градиентного спуска на примере квадратичных функций. Для этого придумайте две-три квадратичные функции от двух переменных, на которых работа методов будет отличаться.

Для анализации траектории градиентного спуска предлагается:

1. Аналитический подход: если функция f(x) является квадратичной, то ее градиент и гессиан также будут квадратичными функциями. Это позволяет аналитически вычислить точное решение градиентного спуска и построить траекторию спуска в виде параболы.

2. Графический подход: можно построить график функции f(x) и нарисовать на нем линии уровня (контуры), соответствующие различным значениям функции. Тогда траектория градиентного спуска будет представлять собой последовательность точек, которые пересекают эти линии уровня.

3. Численный подход: можно запустить градиентный спуск на функции f(x) и записывать значения x на каждой итерации. Затем можно построить график траектории спуска, который будет отображать изменение значения x на каждой итерации. Если функция f(x) является квадратичной, то траектория спуска будет иметь форму параболы.

### Для каждой функции...

#### ... исследуйте сходимость градиентного спуска с постоянным шагом, сравните полученные результаты для выбранных функций

http://www.machinelearning.ru/wiki/index.php?title=%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BD%D0%BE%D0%B3%D0%BE_%D1%81%D0%BF%D1%83%D1%81%D0%BA%D0%B0#:~:text=%D0%A1%D1%85%D0%BE%D0%B4%D0%B8%D0%BC%D0%BE%D1%81%D1%82%D1%8C%20%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BD%D0%BE%D0%B3%D0%BE%20%D1%81%D0%BF%D1%83%D1%81%D0%BA%D0%B0%20%D1%81%20%D0%BF%D0%BE%D1%81%D1%82%D0%BE%D1%8F%D0%BD%D0%BD%D1%8B%D0%BC%20%D1%88%D0%B0%D0%B3%D0%BE%D0%BC,-%D0%A2%D0%B5%D0%BE%D1%80%D0%B5%D0%BC%D0%B0%201%20%D0%BE&text=%D1%80%D0%B5%D0%B0%D0%BB%D0%B8%D0%B7%D1%83%D0%B5%D1%82%D1%81%D1%8F%20%D1%81%D0%B5%D0%B4%D0%BB%D0%BE%2C%20%D0%B0%20%D0%BD%D0%B5%20%D0%BC%D0%B8%D0%BD%D0%B8%D0%BC%D1%83%D0%BC,%D0%BD%D0%B0%D1%85%D0%BE%D0%B4%D1%8F%D1%82%20%D0%BB%D0%BE%D0%BA%D0%B0%D0%BB%D1%8C%D0%BD%D1%8B%D0%B5%20%D0%BC%D0%B8%D0%BD%D0%B8%D0%BC%D1%83%D0%BC%D1%8B%20%D1%86%D0%B5%D0%BB%D0%B5%D0%B2%D0%BE%D0%B9%20%D1%84%D1%83%D0%BD%D0%BA%D1%86%D0%B8%D0%B8.

#### ... сравните эффективность градиентного спуска с использованием одномерного поиска с точки зрения количества вычислений минимизируемой функции и ее градиентов

#### ... исследуйте работу методов в зависимости от выбора начальной точки

#### ... исследуйте влияние нормализации (scaling) на сходимость на примере масштабирования осей плохо обусловленной функции

#### ... в каждом случае нарисуйте графики с линиями уровня и траекториями методов

https://matplotlib.org/

### Реализуйте генератор случайных квадратичных функций n переменных с числом обусловленности k

#### [Сингулярное разложение (попробуйте сначала это объяснить)](https://en.wikipedia.org/wiki/Singular_value_decomposition)

https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html

#### По мнению Chat GPT

Число обусловленности для квадратичной функции определяется как отношение максимального и минимального собственных значений гессиана функции. Для функции f(x,y) = ax^2 + bxy + cy^2, где a, b, c - константы, число обусловленности равно |(a+c+√((a-c)^2 + b^2))/b|. Если b=0, то число обусловленности равно |a/c|. Число обусловленности показывает, насколько чувствительна функция к изменению входных параметров и может быть использовано для оценки устойчивости методов оптимизации. Чем больше число обусловленности, тем менее устойчивым будет метод оптимизации.

### Исследуйте зависимость числа итераций ...

### Для получения более корректных результатов проведите множественный эксперимент и усредните полученные значения числа итераций
