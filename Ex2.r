# this is assigment 4 for the stats class - 01/12/2021
# author: Lorenzo Beltrame

# libs
library(glmnet)
library(selectiveInference)

# set the seed
set.seed(42)

# set working directory
setwd("C:/Users/beltr/Documents/GitHub/Stats")



# preparing the datset
data <- read.csv(file = "communities.csv", header = FALSE)
n <- dim(data)[1]
d <- dim(data)[2]
y <- data[, d]
# remove non-predictor columns (metadata, i.e. community identifier)
X <- data[, -c(1:5, d)]
p <- dim(X)[2]
X <- matrix(as.numeric(as.matrix(X)), nrow = n)
bool <- rowSums(is.na(X)) == 0
X <- X[bool, ]
y <- y[bool]

# adding the names of the columns
colnames <- read.csv(file = "names.csv", header = FALSE, sep = "\n")
# get rid of useless metadata
colnames <- as.character(unlist(colnames))[-(1:5)]
colnames <- sub("@attribute ", "", colnames)
colnames <- sub(" numeric\\\\", "", colnames)
colnames <- sub(" numeric", "", colnames)
colnames <- sub(" numeric}", "", colnames)[1:p]
colnames(X) <- colnames


# I decide to work with less samples
n <- 100
subsample <- sample(1:nrow(X), n)
X <- X[subsample, ]
y <- y[subsample]
p <- dim(X)[2]


# LASSO model selection

lambda <- 0.058
# choose to have only the best 7 features
M <- glmnet(X, y, lambda = lambda)
coeffs <- coef(M)
# print only the non null coeff
chosen_coef <- colnames[coeffs@i]
print("The selected features are: ")
print(chosen_coef) # selected variables

# naive inference
# fit a linear model (OLS) using only the selected variables
fitM <- lm(y ~ X[, chosen_coef])
print(summary(fitM))
# post-model-selection estimate of the error standard deviation (overestimates the true unknown error variance)
sigma <- summary(fitM)$sigma
print(paste0("Sigma is: ", sigma))

# conditional PoSI
lassoInference <- fixedLassoInf(
    x = X,
    y = y,
    family = "gaussian",
    beta = coef(M, s = lambda / n)[-1],
    lambda = lambda,
    alpha = 0.05
)
print(lassoInference)


# uniform PoSI
B.PoSI <- function(alpha, q, N, Nsim = 10000, t.eps = 0.001) {
    G2 <- rchisq(Nsim, df = q)
    t <- seq(0, max(sqrt(G2)), by = t.eps)
    E_G <- 1
    i <- 1
    while (E_G[i] > alpha) {
        mins <- pmin(1, N * (1 - pbeta(t[i]^2 / G2, shape1 = 1 / 2, shape2 = (q - 1) / 2)))
        E_G <- c(E_G, mean(mins))
        i <- i + 1
    }
    # plot(t[1:i], E_G, type="l")
    # abline(h=alpha)
    return(t[i])
}


# if we want to protect against ANY model selection procedure, we should use
kappa <- sum((1:p + 1) * choose(p, 1:p))
result <- B.PoSI(0.05, min(kappa, p, n), kappa, Nsim = 10000, t.eps = 0.01)

print(paste0("The resulting B.PoSI: ", result))


# adjust the naive inference
summary(fitM)
# Compare absolute t-statistics to critical value B.PoSI. Only var45 is still significant. When protecting against ANY model selection procedure, not even var45 is significant any more.
colnames[c(45)] 