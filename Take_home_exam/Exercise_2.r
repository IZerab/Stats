# Lorenzo Beltrame's submission for the 2021WS Stats take home project. 08/02/22

# standard libraries
library(ggplot2)
library(Rmisc)
library(VGAM)
library(jmuOutlier)

set.seed(42)

# my functions
a_f <- function(x, k, m) {
    # If no changes are made
    if (k == 0) {
        return(local_sensitivity(x))
    }
    # Consider the adaption of x
    n <- length(x)

    x <- sort(x)
    x <- c(rep(0, length(x)), x, rep(m, length(x)))

    # Compute point associated to the median
    midpoint <- ceiling(n / 2)

    # Make indices as defined before
    indices <- get_indices(midpoint, k)

    # Vector to store results
    ls <- numeric(length(indices))

    # For each set of indices
    for (i in seq_len(length(indices))) {
        # Extract indices
        ind1 <- indices[[i]][1]
        ind2 <- indices[[i]][2]

        # Compute potential maximum of change of median
        l1 <- x[ind1] - x[ind1 - 1]
        l2 <- x[ind2] - x[ind1]
        l3 <- x[ind2 + 1] - x[ind2]

        # Maximum
        ls[i] <- max(l1, l2, l3)
    }
    # Maximum over all k
    result <- max(ls)
    # Return results
    return(result)
    }

a_fs <- function(x, m) {
    # This function finds the smallest k that follows a_fs > b
    # Number of instances
    n <- length(x)

    # To store computations
    result <- numeric(n + 1)
    names(result) <- paste0("k=", 0:n)

    # Maximum change in median for each k number of changes
    for (k in 0:n) {
        result[k + 1] <- a_f(x, k, m)
    }

    # Return results
    return(result)
    }

d_f <- function(a, b) {
    # function that actually computes d_f
    result <- which(a > b)[1]
    return(result - 1)
    }

propose_test_release <- function(x, a, b, m, size=1000) {
    # Initialize result
    result <- numeric(size)

    # For each iteration
    w <- rlaplace(n = size) # nolint
    d <- d_f(a, b)

    # Left-hand side and right-hand side
    left <- d + w / epsilon
    right <- log(2 / delta) / (2 * epsilon)
    output <- median(x) + (b / epsilon) * rlaplace(size) # nolint
    result <- ifelse(left < right, NA, output)
    return(result)
    }

local_sensitivity <- function(x) {
    # Function that computes the local sensitivity of the data passed
    # sort x in order
    x <- sort(x)

    # Number of instances considered
    n <- length(x)

    # compute the point of the median
    midpoint <- ceiling(n / 2)

    # Compute the local sensitivity
    temp <- x[midpoint + 1] - x[midpoint]
    temp1 <- x[midpoint] - x[midpoint - 1]
    result <- max(temp, temp1)

    return(result)
    }

inverse_local_sensitivity <- function(x, m) {
    # This function computes the inverse local sensitivity for the
    # data passed in x

    # First I sort x
    x <- sort(x)
    n <- length(x)

    # Since n is odd, the middle point is:
    mid <- ceiling(n / 2)

    # Initialise the output
    result <- numeric(m + 1)
    names(result) <- paste0("z=", 0:m)

    # For each potential value of the median
    for (i in 1:(m + 1)) {
        temp <- x
        # Potential median
        z <- i - 1

        # While the current median is not equal to z
        while (median(temp) != z) {
            # z is the new middle point
            temp[mid] <- z

            # Sort the data again
            temp <- sort(temp)

            # Add count
            result[i] <- result[i] + 1
        }
    }

    return(result)
    
    }

get_indices <- function(mid, max_entries_changed) {
    # This function gets the list of potential mid points
    results <- list()

    # support variable
    counter <- 1

    # For each number of potential changes k
    for (i in 1:max_entries_changed) {
        # Compute the number of way that the changes might influence the median
        for (j in 1:i) {
            results[[counter]] <- c(mid - i, mid + 1) + (j - 1)
            counter <- counter + 1
        }
    }
    # Return results
    return(results)
}

pmf_z <- function(inv_local, epsilon) {
    # Function that computes the probability mass function of Z based
    # on the inverse local sensitivity

    num <- exp((-epsilon / 2) * inv_local)
    den <- sum(exp((-epsilon / 2) * inv_local))

    # Compute the fraction
    result <- num / den
    # Return results
    return(result)
}

# Generate the required data and model hyperparameters

# DP hyperparameters
epsilon <- 1
delta <- 0.01
b <- 10

# data
n <- 101
max_value <- 50
x <- sample(0:max_value, size = n, replace = TRUE, prob = 1:(max_value + 1))

# plot the data sampled
dev.new()
hist(x)

# True median
cat(paste0("True median = ", median(x), "\n"))

# Compute the inverse local sensitivity
inverse_local <- inverse_local_sensitivity(x, max_value)
print(inverse_local)

# Compute pmf of Z
q <- pmf_z(inverse_local, epsilon)

# Plot pmf and actual median
dev.new()
pmf_plot <- plot(0:max_value, q, type = "b", pch = 4, xlab = "Z", ylab = "prob")
abline(v = median(x), lty = 2, col = "#2c8d3c")
text(paste0("True median = ", median(x)), x = 3, y = max(q), adj=c(0,1))


# Sample 1000 values according to pmf of Z
size <- 1000
median_to_publish <- sample(0:max_value, size = size, replace = TRUE, prob = q)

# For later
med_ils <- median_to_publish

# Plot histogram of published medians
breaks <- seq(min(median_to_publish), max(median_to_publish), length.out = 10)
dev.new()
h <- hist(median_to_publish, xlim = c(0, 50), breaks = breaks,
    main = "Inverse local sensitivity", xlab = "Median to publish")

text(paste0("True median = ", median(x)), x = 0.1, y = max(h$counts),adj = c(-0.5, 1))


# My proposal for a test release
n <- 101
m <- 50
epsilon <- 1
delta <- 0.01
prop_test_relase <- list()

x <- sample(0:m, size = n, replace = TRUE, prob = 1:(m + 1))
# True median
cat(paste0("True median = ", median(x)))
cat("\n")

# Try different values for b

list_b <- c(0.5, 0.7, 1, 2, 3, 4, 5, 6, 7, 8)
propose_test_relase <- list()
a <- a_fs(x, m)

# Compute propose-test-release
for (i in seq_len(length(list_b))) {
    b <- list_b[i]
    prop_test_relase[[i]] <- propose_test_release(x, a, b, m)
}

# print the percentage of NaN in the proposed test released
cat("Percentage of NaNs: \n")
percentage <- 100 * sapply(lapply(prop_test_relase, is.na), mean)
print_temp <- paste0("b=", list_b, ": ", percentage, "% are NaNs")

for (t in print_temp) {
      print(t)
}

# Plot the results
xmin <- min(0, sapply(prop_test_relase, min, na.rm = TRUE))
xmax <- max(m, sapply(prop_test_relase, max, na.rm = TRUE))
breaks <- seq(xmin, xmax, length.out = 80)
maxcounts <- numeric(length(list_b))
for (i in seq_len(length(list_b))) {
    b <- list_b[i]
    p <- prop_test_relase[[i]]
    h <- hist(p, breaks = breaks, plot = FALSE)
    maxcounts[i] <- max(h$count)
    }
    ymax <- max(maxcounts)
    par(mfrow = c(3, 1), mar = c(2, 2, 1, 1))
    for (i in seq_len(length(list_b))) {
    b <- list_b[i]
    p <- prop_test_relase[[i]]
    
    
    dev.new()
    hist(
        p,
        xlim = c(0, 100),
        ylim = c(0, ymax),
        breaks = breaks,
        xlab = "Median to publish",
        main = ""
        )
    
    abline(v=median(x), col="blue")
    
    text(labels = paste0("True median = ", median(x)), x = 3, y = ymax, adj=c(-1.5,1))
    
    mtext(paste0("b=", list_b[i]) )
    
    
}

#par(mfrow = c(3, 1), mar = c(2, 2, 1, 1))
