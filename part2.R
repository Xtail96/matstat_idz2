library(maxLik)
# шаг, заданный в условии
m_h <- 1.10

# число элементов в выборке
m_n <- 50

# левая граница промежутка попадания
m_c <- 1.43

# правая граница промежутка попадания
m_d <- 3.13

m_alpha2 <- 0.05
m_lambda0 <- 0.36
m_lambda1 <- 0.21

# функция для нахождения эмпирической функции
empirical.function <- function(x, t)
{
  z <- x[x < t]
  length(z) / length(x)
}

# ------------------------------ begin ----------------------------------------
print("Part2")

# моя выборка из задания
m_selection <- scan("input_part2.txt")
print("Исходные данные:")
print(m_selection)
print("-------------------------------------------------------")

# ----- a -----

# вариационный ряд
selection.sorted <- sort(m_selection)
print("Вариационный ряд:")
print(selection.sorted)
print("-------------------------------------------------------")

# содержит в себе уникальные аргументы выборки
selection.unique.x <- unique(selection.sorted)

# значения эмпирической функции распределения в каждом уникальном аргументе выборки
selection.unique.y <- 0

for(i in 1:length(selection.unique.x))
{
  selection.unique.y[i] <- empirical.function(m_selection, selection.unique.x[i])
  selection.unique.y[length(selection.unique.x) + 1] <- 1
}

# построение ступенчатого вида эмпирической функции распределения
selection.empirical.function <- stepfun(selection.unique.x, selection.unique.y)
plot(selection.empirical.function, col="blue", verticals = TRUE, main = "Эмпирическая функция распределения:")

# верхняя граница (целые числа, не меньшие соответсвующих чисел в переданом векторе, в данном случае вектор состоит из 1 элемента)
# в качестве верхней границы берем верхнюю границу группы с максимальным в выборке значением
upper.bound <- ceiling(selection.unique.x[length(selection.unique.x)])
#upper.bound <- ceiling(selection.unique.y[m_h] - selection.unique.y[1]/m_h)

# вектор границ для гисограммы
# заполняем первый элемент вручную
k <- selection.sorted[1] - m_h/2
for(i in 2:(upper.bound + 2))
{
  k[i] <- k[i-1] + m_h
}

h <- hist(selection,
          col="lightblue",
          main = "Гистограмма частот",
          xlab = "Элементы выборки",
          ylab = "Частота встречи",
          breaks = k,
          right = TRUE)

# ----- b -----
print("Выборочное среднее (Математическое ожидание):")
selection.avg <- sum(m_selection) / length(m_selection)
print(selection.avg)
print("-------------------------------------------------------")

print("Выборочная Дисперсия:")
s2 <- (1/length(m_selection)*sum((m_selection - selection.avg)^2))
print(s2)
print("-------------------------------------------------------")

print("Выборочная медиана:")
median <- sort(m_selection)[trunc(length(m_selection)*1/2+1)]
print(median)
print("-------------------------------------------------------")

print("Выборочная ассиметрия:")
assymetry <- (1/length(m_selection)) * sum((m_selection - selection.avg)^3) / s2^(3/2)
print(assymetry)
print("-------------------------------------------------------")

print("Выборочный эксцесс:")
exc <- ((1/length(m_selection)) * sum((m_selection - selection.avg)^4) / s2^2) - 3
print(exc)
print("-------------------------------------------------------")

print("Вероятность попадания в промежуток [с; d]:")
print("c:")
print(m_c)
print("d:")
print(m_d)
p <- empirical.function(m_selection, m_c) - empirical.function(m_selection, m_d)
print("p:")
print(p)
print("-------------------------------------------------------")


