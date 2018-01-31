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
upper.bound <- ceiling(max(m_selection))

# вектор границ для гисограммы
# заполняем первый элемент вручную
k <- 0;
for(i in 2:(upper.bound))
{
  k[i] <- k[i-1] + m_h
}

h <- hist(m_selection,
          col="lightblue",
          main = "Гистограмма частот",
          xlab = "x",
          ylab = "Частота встречи",
          breaks = k,
          right = FALSE)

plot(h$counts ~ h$mids,
     col="red",
     type="l",
     bty="n",
     main="Полигон частот для x",
     xlab="x",
     ylab="Частота")

# ----- b -----
print("Выборочное среднее (Математическое ожидание):")
selection.mean <- sum(m_selection) / length(m_selection)
print(selection.mean)
print("-------------------------------------------------------")

print("Несмещенная Выборочная Дисперсия:")
s2 <- (1/(length(m_selection) - 1)*sum((m_selection - selection.mean)^2))
print(s2)
print("-------------------------------------------------------")

print("Выборочная медиана:")
median <- sort(m_selection)[trunc(length(m_selection)*1/2+1)]
print(median)
print("-------------------------------------------------------")

print("Выборочная ассиметрия:")
assymetry <- (1/length(m_selection)) * sum((m_selection - selection.mean)^3) / s2^(3/2)
print(assymetry)
print("-------------------------------------------------------")

print("Выборочный эксцесс:")
exc <- ((1/length(m_selection)) * sum((m_selection - selection.mean)^4) / s2^2) - 3
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

# ----- c -----
print("Оценка максимального правдоподобия:")
log.lik <- function(theta)
{
  sum(log(theta) - theta * m_selection)
}
grad.lik <- function(theta)
{
  sum(1/theta - m_selection)
}
hess.lik <- function(theta)
{
  -100/theta^2
}
a <- maxNR(log.lik, grad.lik, hess.lik, start = 1)
#print(summary(a))
print("a$estimate:")
print(a$estimate)
print("a$gradient:")
print(a$gradient)
print("-------------------------------------------------------")

# ----- d -----
print("Асимптотический доверительный интервал, уровня значимости alpha2")
# ОМП +- (x.alpha/sqrt(n*I(ОМП)))
print("alpha2:")
print(m_alpha2)
print("Интервал T:")
omp <- a$estimate
x.alpha <- qnorm(1 - m_alpha2/2)
I = omp^2
d <- x.alpha / sqrt(length(m_selection) * I)
T <- array(dim = 2)
T[1] <- selection.mean - d
T[2] <- selection.mean + d
print(T)
print("-------------------------------------------------------")

# ----- e -----
print("Критерий значимости проверки простой гипотезы с использованием теоремы Колмогорова:")
#a0 <- selection.mean
#s0 <- s2
#x <- m_selection
#n <- length(m_selection)

# скрипт из методички с изменением распределения
#v1 <- sort(x)
#v2 <- c(0:(n-1))/n
#v3 <- c(1:n)/n
#v4 <- abs(pexp(v1, a0, s0) - v2)
#v5 <- abs(pexp(v1, a0, s0) - v3)
#D <- max(v4, v5)
#print("D:")
#print(D)

# проверям простую гипотезу H0 согласия с показательным распределением с параметрами lambda0
print("lambda0:")
print(m_lambda0)

# генерируем показательное распределение
#kolm <- rexp(n, m_lambda0)
# проверяем нашу выборку на соответствие сгенерированному выше распределению
#test1 <- ks.test(m_selection, kolm)
#print(test1)

# проверяем нашу выборку на соответствие показательному распределению без явной генерации тестового распределения
test2 <- ks.test(m_selection, "pexp", m_lambda0)
print(test2)

max.value <- test2$p.value
print("max_value (p-value):")
print(max.value)
print("-------------------------------------------------------")

# ----- f -----
# выделяем промежутки ([0; 1.1) [1.1; 2.2)  [2.2; 4.4) [4.4; 8.8) [8.8; +inf) )
# считаем практическое значение
print("nk (значения, которые получились):")
nk <- c(25, 14, 3+3, 0+1+0+2, 0+1+0+0+1)
r <- 5
print(nk)

# считаем теоретическое значение
print("pk (плотность каждого интервала):")
pk <- array(dim = 5)
left.bound <- c(0.0, 1.1, 2.2, 4.4)
right.bound <- c(1.1, 2.2, 4.4, 8.8)
pk <- pexp(right.bound, m_lambda0) - pexp(left.bound, m_lambda0)
pk[5] <- 1 - sum(pk)
print(pk)

print("npk (значения, которые должны были получиться):")
npk <- pk*n
print(npk)

print("test1:")
test1 <- chisq.test(nk, p = pk)
print(test1)

print("xi^2 экспериментальное:")
xi2 <- as.numeric(test1[1]$statistic)
print(xi2)

print("Теоретическое значение на уровне alpha2:")
print("alpha2:")
print(m_alpha2)
print("value:")
kvantil <- qchisq(1 - m_alpha2, 4)
print(kvantil)

if(xi2 > kvantil)
{
  print("Гипотезу H0 нужно отвергнуть")
}else
{
  print("Гипотезу H0 нужно принять")
}

print("maxValue:")
max.value <- 1 - pchisq(xi2, 4)
print(max.value)

print("-------------------------------------------------------")

# ----- g -----
# функция, минимум которой нужно найти
theta <- function(lambda)
{
  pk <- pexp(right.bound, lambda) - pexp(left.bound, lambda)
  pk[5] <- 1 - sum(pk)
  sum(((length(m_selection)*pk - nk)^2)/(length(m_selection)*pk))
}
xm <- nlm(theta, p = mean(selection.sorted))
min.lambda <- xm$estimate
print("Минимальное значение lambda:")
print(min.lambda)


print("pk (что должно было получиться):")
pk <- pexp(right.bound, min.lambda) - pexp(left.bound, min.lambda)
pk[5] <- 1 - sum(pk)
print(pk)

print("test1:")
test1 <- chisq.test(nk, p = pk)
print(test1)

print("xi2 экспериментальное:")
xi2 <- as.numeric(test1[1]$statistic)
print(xi2)


print("Теоретическое значение на уровне alpha2:")
print("alha2:")
print(m_alpha2)
print("value:")
kvantil <- qchisq(1 - m_alpha2, 3)
print(kvantil)

if(xi2 > kvantil)
{
  print("Гипотезу H0 нужно отвергнуть")
}else
{
  print("Гипотезу H0 нужно принять")
}

print("maxValue:")
max.value <- 1 - pchisq(xi2, 3)
print(max.value)
print("-------------------------------------------------------")

# ----- h -----
prod0 = prod(dexp(m_selection, m_lambda0))
prodA = prod(dexp(m_selection, m_lambda1))

print("test1. Для основной гипотезы:")
print("Отношение правдоподобия:")
ratio <- prod0/prodA
print(ratio)
print("Критическое значение:")
critical.value <- qexp(1 - m_alpha2, m_lambda1)
print(critical.value)
if(ratio > critical.value)
{
  print("отклоняем основную гипотезу (X ~ Exp(lambda)) и принимаем альтернативу")
} else
{
  print("принимаем гипотезу X ~ Exp(lambda) на уровне значимости alpha2")
}

print("Меняем гипотезы местами")

print("test2. Для для альтернативной гипотезы:")
print("Отношение правдоподобия:")
ratio <- prodA/prod0
print(ratio)
print("Критическое значение:")
critical.value <- qexp(1 - m_alpha2, m_lambda0)
print(critical.value)

if(ratio > critical.value)
{
  print("отклоняем основную гипотезу (X ~ Exp(lambda)) и принимаем альтернативу")
} else
{
  print("принимаем гипотезу X ~ Exp(lambda) на уровне значимости alpha2")
}
print("-------------------------------------------------------")

