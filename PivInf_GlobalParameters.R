################################################
# Global parameters throughout all simulations #
################################################

sigma        <- 1                      # sd of innovations
n_disc       <- 20                     # Number discretization points of [0,1]
lam_grid     <- (1:n_disc)/n_disc      # Values of lambda for seq. statistics
trunk        <- 200                    # Lin. process approximated by MA(trunk) 
Delta        <- seq(0, 1, by = 0.01)   # Thresholds for emp. rej. probs
sample_sizes <- c(100, 200, 500, 1000) # Sample sizes
rep          <- 1000                   # Number of repetitions