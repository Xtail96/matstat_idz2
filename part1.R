library(maxLik)

selection <- scan("input_part1.txt")
print("Исходные данные:")
print(selection)

# a
s <- c(0:8)
h <- hist(selection, col="lightgreen",
     main = "Гистограмма частот",
     xlab = "Элементы выборки",
     ylab = "Частота встречи",
     breaks = s,
     right = F)

selection.sorted <- sort(selection)
print("Вариационный ряд:")
print(selection.sorted)

f <- function(x, t)
{
  z <- x [x < t]
  length(z)/length(x)
}

xu <- unique(selection.sorted)
yu <- 0

for(i in 1:length(xu))
{
  yu[i] <- f(selection, xu[i])
  yu[length(xu) + 1] <- 1
}

z <- stepfun(xu, yu)
plot(z, col="blue", verticals = TRUE, main = "Эмпирическая функция распределения:")

#v <- seq(from <- 1/50, to <- 1, by <- 1/50)
#plot(selection.sorted, v)

# b
mean <- mean(sum(selection)/length(selection))
print("Выборочное среднее:")
print(mean)

disp <- sum(selection^2) / length(selection) - mean^2
print("Выборочная дисперсия:")
print(disp)

median <- selection.sorted[trunc(length(selection) * 1/2 + 1)]
print("Выборочная медиана:")
print(median)

asymmetry <- 1/length(selection) * ( sum((selection - mean)^3) / disp^(3/2) )
print("Выборочная ассиметрия:")
print(asymmetry)

excess <- 1/length(selection) * ( sum((selection - mean)^4) / disp^2 ) - 3
print("Выборочный эксцесс:")
print(excess)

a <- 0.00
b <- 2.41
propabitity.of.hitting <- f(selection, b) - f(selection, a)
print("Вероятность попадания в промежуток [a = 0.00, b = 2.41]:")
print(propabitity.of.hitting)

# c
LL <- function(t)
{
  sum(dpois(selection, t[1], log = TRUE))
}

ml <- maxNR(LL, start = c(1))
val <- ml$estimate
print("Смещение оценок:")
print(val)

# d
mle<-val;
a1 <- 0.05
xa1 <- qnorm(1-a1/2)
T <- array(dim = 2)
T[1] <- mle - xa1 * sqrt(mle/length(selection))
T[2] <- mle + xa1 * sqrt(mle/length(selection))
print("Асимптотический доверительный интервал уровня значимости a1 = 0.05:")
print(T)

#e
nk <- c(20, 19, 7, 2, 2)
print("nk:")
print(nk)
pk <- dpois(c(0, 1, 2, 3), 2)
print("pk:")
pk <- c(pk, (1 - sum(pk)))
print(pk)
xi2 <- as.numeric(chisq.test(nk, p = pk)[1]$statistic)
kvantil <- qchisq(1 - a1, 4)

print("xi2:")
print(xi2)
print("kvantil:")
print(kvantil)

if(xi2 > kvantil)
{
  print("Гипотезу H0 нужно отвергнуть")
}else
{
  print("Гипотезу H0 нужно принять")
}

max.value <- 1 - pchisq(xi2, 4)
print("maxValue:")
print(max.value)

#f
teta <- function(lambda)
{
  pk <- dpois(c(0, 1, 2, 3), lambda)
  pk <- c(pk, (1 - sum(pk)))
  sum(((length(selection)*pk - nk)^2)/(length(selection)*pk))
}
xm <- nlm(teta, p = mean(selection.sorted))
min.lambda <- xm$estimate
print("minLambda:")
print(min.lambda)


pk <- dpois(c(0, 1, 2, 3), min.lambda)
pk <- c(pk, (1 - sum(pk)))

print("pk:")
print(pk)

xi2 <- as.numeric(chisq.test(nk, p = pk)[1]$statistic)
kvantil <- qchisq(1 - a1, 3)

print("xi2:")
print(xi2)
print("kvantil:")
print(kvantil)

if(xi2 > kvantil)
{
  print("Гипотезу H0 нужно отвергнуть")
}else
{
  print("Гипотезу H0 нужно принять")
}

max.value <- 1 - pchisq(xi2, 3)
print("maxValue:")
print(max.value)

#csq <- function (t)
#{
#  p <- ppois(k[2:(r)], t) - ppois(k[1:r-1], t)
  #print(nu)
  #print(p)
#  f <- sum((nu-n*p)^2/(n*p))
  #print(p)
  #print(f)
#}

#xi2 <- nlm(csq, p = mean(selection))
#xia11 <- qchisq(1 - a1, r-2)

#print("xi^2:")
#print(xi2$minimum)
#print("xi alpha 11:")
#print(xia11)
#if(xi2$minimum <= xia11)
#{
#  print("Гипотезу H0 нужно принять")
#}else
#{
#  print("Гипотезу H0 нужно отвергнуть")
#}
