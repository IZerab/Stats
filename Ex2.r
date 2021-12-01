# this is assigment 4 for the stats class - 01/12/2021
# author: Lorenzo Beltrame

# libs
library(glmnet)
library(selectiveInference)

# set working directory
setwd("/Users/Lorenzo/GitHub/Stats")

# preparing the datset
data <- read.csv("communities.data", header = FALSE)
n <- dim(data)[1]
p <- dim(data)[2]
Y <- data[,d]

# in my design matrix I do not wat community identifiers!
X <- data[, -c(1:5, d)]

p <- dime(X)[2]
X <- matrix(as.numeric(as.matrix(X)), nrow = n)
bool <- rowSums(is.na(X)) == 0
X <- X[bool,]
Y <- Y[bool,]


# get the names of the columns (have been saven in a separated file!)
colnames <- read.table(file="names.csv", header=FALSE, sep="\n")
colnames <- as.character(unlist(colnames))[-(1:11)]

# clean what I imported via substitution
colnames <- sub("@attribute ", "", colnames)
colnames <- sub(" numeric\\\\", "", colnames)
colnames <- sub(" numeric}", "", colnames)[1:p]

colnames(X) <- colnames


# I decide to work with less samples
n <- 100
set.seed(42)

subsample <- sample(1:nrow(X), n)
X <- X[subsample]
Y <- Y[subsample]
p <- dim(X)[2]





graphics.off()  # clear all graphs
    rm(list = ls()) # remove all files from your workspace
    
    n <- 500 # number of observations
    p <- 20  # number of variables


# FIRST TASK: naive approach

