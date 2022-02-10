# Lorenzo Beltrame's submission for the 2021WS Stats take home project. 08/02/22

# standard libraries
library(ggplot2)
library(Rmisc)



# number of samples
compute_power <- function(alpha, n=1000, mmax=6) {
    # X-Axis spacing (the fixed value of the statistic)
    xx <- seq(0, mmax, length.out = n)

    # Compute the power for A and Corrected p-value
    pp_a <- numeric(n)
    pp_cor <- numeric(n)

    for (i in 1:n) {
        # Copute power of p_a
        arg <- qnorm(1 - alpha) - xx[i]
        pp_a[i] <- 1 - pnorm(arg)
        # Compute power of p_cor
        arg1 <- qnorm(1 - alpha / 2) - xx[i]
        arg2 <- - qnorm(1 - alpha / 2) - xx[i]
        pp_cor[i] <- 1 - pnorm(arg1) + pnorm(arg2)
    }

    # Return the results
    res <- data.frame(xx = xx, pp_a = pp_a, pp_cor = pp_cor)
    return(res)
}


# run the experiments!

# values of alpha chosen
alphas <- c(0.1, 0.05, 0.01, 0.001, 0.000001, 0.000000001)
values <- list()
plots <- list()


for (i in seq_len(length(alphas))) {
    values[[i]] <- compute_power(alphas[i])
}

# Plotting

for (i in seq_len(length(alphas))) {

    # convenience variables
    xx <- values[[i]]$xx
    pp_a <- values[[i]]$pp_a
    pp_cor <- values[[i]]$pp_cor

    # dataframe for the plots
    data <- data.frame(index = xx, a_power = pp_a, cor_power = pp_cor)

    # create the plot
    string <- "Corrected"
    title <- paste("Power of the ", as.character(alphas[i]), " test")
    my_plot <- ggplot(data, aes(x = index)) +
                geom_line(aes(y = a_power, colour = "'A' test"), size = 0.9) +
                geom_line(aes(y = cor_power, colour = string), size = 0.9) +
                labs(colour = "Type of test") +
                ylim(0, 1) +
                theme(
                    legend.position = c(1, .4),
                    legend.justification = c("right", "top"),
                    legend.box.just = "right",
                    legend.margin = margin(2, 2, 2, 2),
                    legend.background = element_rect(fill = "#eee8e898")) +
                labs(title = title,
                     x = "Values of Delta", y = "Power") +
                theme(plot.title = element_text(hjust = 0.5))

    plots <- append(plots, list(my_plot))
}

multiplot(plotlist = plots, cols = 2)
