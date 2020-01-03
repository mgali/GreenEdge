# Test CI of pre-screening test using bootstrapping

data <- read.csv(file = "~/Desktop/GreenEdge/GCMS/pre-screening_test.csv")

summary(data$ratio)

nboot <- 100000

# Untransformed data
a <- numeric(nboot)
for (i in 1:nboot) {
  a[i] <- mean(sample(data$ratio, replace = T))
}

summary(a)
hist(a, plot = T)

quantile(a, c(0.025, 0.975))
quantile(a, c(0.005, 0.995))

# Log ratio
b <- numeric(nboot)
for (i in 1:nboot) {
  b[i] <- mean(sample(log10(data$ratio), replace = T))
}

summary(b)
hist(b, plot = T)

quantile(b, c(0.025, 0.975))
quantile(b, c(0.005, 0.995))