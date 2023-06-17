# Лабораторная работа №1. Градиентный спуск

## Создатели

* [@sfbakturin](https://github.com/sfbakturin)
* [@Nomad192](https://github.com/Nomad192)
* [@SotnikovMaksim](https://github.com/SotnikovMaksim)

## Условие

### Постановка задачи

1. Реализуйте градиентный спуск с постоянным шагом (learning rate).
2. Реализуйте метод одномерного поиска (метод дихотомии, метод Фибоначчи, метод золотого сечения) и градиентный спуск на его основе.
3. Проанализируйте траекторию градиентного спуска на примере квадратичных функций. Для этого придумайте две-три квадратичные функции от двух переменных, на которых работа методов будет отличаться.
4. Для каждой функции:
    * исследуйте сходимость градиентного спуска с постоянным шагом, сравните полученные результаты для выбранных функций;
    * сравните эффективность градиентного спуска с использованием одномерного поиска с точки зрения количества вычислений минимизируемой функции и ее градиентов;
    * исследуйте работу методов в зависимости от выбора начальной точки;
    * исследуйте влияние нормализации (scaling) на сходимость на примере масштабирования осей плохо обусловленной функции;
    * в каждом случае нарисуйте графики с линиями уровня и траекториями методов;
5. Реализуйте генератор случайных квадратичных функций *n* переменных с числом обусловленности *k*.
6. Исследуйте зависимость числа итераций *T(n, k)*, необходимых градиентному спуску для сходимости в зависимости от размерности пространства 2 ⩽ *n* ⩽ 10^3 и числа обусловленности оптимизируемой функции 1 ⩽ *k* ⩽ 10^3.
7. Для получения более корректных результатов проведите множественный эксперимент и усредните полученные значения числа итераций.

### Дополнительное задание

Реализуйте одномерный поиск с учетом условий Вольфе и исследуйте его эффективность. Сравните полученные результаты с реализованными ранее методами.

## Структура проекта

### Каталоги проекта

* [**`Image/`**](Image/) - директория с картинками для отчета. Все названия соответствуют следующему шаблону: *Ta_Fc_\*.png* или же *Ta_Pb_Fc_\*.png*, где *a* - номер задания, *b* - номер пункта, *\** - дополнительная информация.
* [**`HQ/`**](HQ/) - директория с высококачественными картинками для просмотра отдельно от отчета. Все названия соответствуют следующему шаблону: *Ta_Pb_Fc_\*_HQ.png*, где *a* - номер задания, *b* - номер пункта, *\** - дополнительная информация.
* [**`Data/`**](Data/) - директория с `.csv` данными от полученных графиков. Все названия соответствуют следующему шаблону: *Ta_Fc_\*.csv* или же *Ta_Pb_Fc_\*.csv*, где *a* - номер задания, *b* - номер пункта, *\** - дополнительная информация.

### Главные файлы в корневой директории

* [**`T1.ipynb`**](T1.ipynb), [**`T2.ipynb`**](T2.ipynb), [**`T3.ipynb`**](T3.ipynb), [**`T4.ipynb`**](T4.ipynb), [**`T5.ipynb`**](T5.ipynb), [**`T6-7.ipynb`**](T6-7.ipynb) - сырые Jupyter Notebook основных заданий.
* [**`AT.ipynb`**](AT.ipynb) - сырой Jupyter Notebook дополнительного задания.
* [**`Report.tex`**](Report.tex) - исходный файл отчета.
* [**`Report.pdf`**](Report.pdf) - собранный файл отчета.
* [**`.gitignore`**](.gitignore) - все игнорируемые `git` файлы.
* [**`README.md`**](README.md) - этот файл.
