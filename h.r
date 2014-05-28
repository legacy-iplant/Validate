# Get data ready
test.data <- read.table(file="/users/dustin/validate/outputplink/plinkstd1.qassoc", header=TRUE)
s <- test.data$P
snps <- test.data$SNP
y <- c(TRUE, TRUE, TRUE, rep(FALSE, length(snps)-3))

# Functions
n <- length(s)
n1 <- sum(y)
n0 <- n-n1
pi0 <- n0/n
pi1 <- n1/n
severity.ratio <- pi0/pi1
zord <- order(s)
sc <- s[zord]

Get.Score.Distributions <- function(y, s, n1, n0) {
  s1 <- unname(tapply(y, s, sum))/n1
  s1 <- c(0, s1, 1 - sum(s1))
  s0 <- unname(tapply(1 - y, s, sum))/n0
  s0 <- c(0, s0, 1 - sum(s0))
  S <- length(s1)
  F1 <- cumsum(s1)
  F0 <- cumsum(s0)
  return(list(F1 = F1, F0 = F0, s1 = s1, s0 = s0, S = S))
}

out.scores <- Get.Score.Distributions(y = y, s = s, n1 = n1, n0 = n0)

F1 <- out.scores$F1
F0 <- out.scores$F0
s0 <- out.scores$s0
s1 <- out.scores$s1
S <- out.scores$S

chull.points <- chull(1 - F0, pmax(1 - F1, 1 - F0))

hc <- length(chull.points)

cost <- c(1:(hc + 1))
b0 <- c(1:hc + 1)
b1 <- c(1:hc + 1)
G0 <- 1 - F0[chull.points]
G1 <- 1 - F1[chull.points]

if (severity.ratio > 0) {
  shape1 <- 2
  shape2 <- 1 + (shape1 - 1) * 1/severity.ratio
}

if (severity.ratio < 0) {
  shape1 <- pi1 + 1
  shape2 <- pi0 + 1
}

cost[1] <- 0
cost[hc + 1] <- 1
b00 <- beta(shape1, shape2)
b10 <- beta(1 + shape1, shape2)
b01 <- beta(shape1, 1 + shape2)
b0[1] <- pbeta(cost[1], shape1 = (1 + shape1), shape2 = shape2) * 
  b10/b00
b1[1] <- pbeta(cost[1], shape1 = shape1, shape2 = (1 + 
                                                     shape2)) * b01/b00
b0[hc + 1] <- pbeta(cost[hc + 1], shape1 = (1 + shape1), 
                    shape2 = shape2) * b10/b00
b1[hc + 1] <- pbeta(cost[hc + 1], shape1 = shape1, shape2 = (1 + 
                                                               shape2)) * b01/b00
for (i in 2:hc) {
  cost[i] <- pi1 * (G1[i] - G1[i - 1])/(pi0 * (G0[i] - 
                                                 G0[i - 1]) + pi1 * (G1[i] - G1[i - 1]))
  b0[i] <- pbeta(cost[i], shape1 = (1 + shape1), shape2 = shape2) * 
    b10/b00
  b1[i] <- pbeta(cost[i], shape1 = shape1, shape2 = (1 + 
                                                       shape2)) * b01/b00
}

LHshape1 <- 0
for (i in 1:hc) {
  LHshape1 <- LHshape1 + pi0 * (1 - G0[i]) * (b0[(i + 
                                                    1)] - b0[i]) + pi1 * G1[i] * (b1[(i + 1)] - b1[i])
}

B0 <- pbeta(pi1, shape1 = (1 + shape1), shape2 = shape2) * 
  b10/b00
B1 <- pbeta(1, shape1 = shape1, shape2 = (1 + shape2)) * 
  b01/b00 - pbeta(pi1, shape1 = shape1, shape2 = (1 + 
                                                    shape2)) * b01/b00
H <- 1 - LHshape1/(pi0 * B0 + pi1 * B1)

