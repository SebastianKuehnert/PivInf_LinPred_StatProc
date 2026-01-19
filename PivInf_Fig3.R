#################################################
# Empirical rejection probabilities in Figure 3 #
#################################################

library(MASS)
library(MTS)

source("PivInf_AuxFunctions.R")
source("PivInf_GlobalParameters.R")

set.seed(123)

# Required parameters for Figure 3 #############################################
d            <- 5
Sigma        <- diag(d)
q            <- 3 # VAR order
alpha        <- 0.1
p_value      <- 1

AR1 <- 0.16 * matrix(c(7, 2, 1, 0, 0,
                       2, 5, 2, 1, 0,
                       1, 2, 5, 2, 1,
                       0, 1, 2, 5, 2,
                       0, 0, 1, 2, 5
              ), nrow = 5, byrow = TRUE)

AR2 <- -0.1 * matrix(c(3, 2, 0, 0, 0,
                       2, 3, 2, 0, 0,
                       0, 2, 3, 2, 0,
                       0, 0, 2, 3, 2,
                       0, 0, 0, 2, 3
              ), nrow = 5, byrow = TRUE)

AR3 <- -0.05 * matrix(c(2, 1, 0, 0, 0,
                        1, 1, 1, 0, 0,
                        0, 1, 1, 1, 0,
                        0, 0, 1, 1, 1,
                        0, 0, 0, 1, 1
               ), nrow = 5, byrow = TRUE)

# Check stationarity of VAR(3) process #########################################
d <- nrow(AR1)
comp_mat <- rbind(
   cbind(AR1, AR2, AR3),
   cbind(diag(d), matrix(0, d, 2*d)),
   cbind(matrix(0, d, d), diag(d), matrix(0, d, d))
)
if (all(Mod(eigen(comp_mat)$values) < 1)) 
   "The VAR(3) process is stationary"

# Coefficients of the linear process
Theta_list <- vector("list", length = trunk + p_value + 1)
Theta_list[[1]] <- diag(d)
for (j in 1:(trunk + p_value)) {
   Theta_list[[j + 1]] <- (0.6 * AR1) %*% Theta_list[[j]]
}

# Covariance matrices ##########################################################

### For VAR(3) process
covs_VAR <- VARMAcov(Phi = cbind(AR1, AR2, AR3), Sigma = Sigma, lag = p_value)[[1]]
Gamma_list_VAR <- lapply(0:(ncol(covs_VAR) / nrow(covs_VAR) - 1), function(i) {
   covs_VAR[, (i * nrow(covs_VAR) + 1):((i + 1) * nrow(covs_VAR))]
})

### For linear vector-valued process
Gamma_list_LIN <- vector("list", length = p_value + 1)
for (h in 0:p_value) {
   Gamma_h <- matrix(0, d, d)
   for (j in 0:trunk) {
      Gamma_h <- Gamma_h + Theta_list[[j + h + 1]] %*% t(Theta_list[[j + 1]])
   }
   Gamma_list_LIN[[h + 1]] <- Gamma_h
}

build_Gpj <- function(p, j, lst) {
   d <- nrow(lst[[1]])
   e <- rep(0, d); e[j] <- 1
   tr_Gamma0 <- as.numeric(t(e) %*% lst[[1]] %*% e)
   G_pj <- matrix(0, d * (p + 1) + 1, d * (p + 1) + 1)
   
   G_pj[1, 1] <- tr_Gamma0
   
   row_blocks <- lapply(1:(p + 1), function(k) t(e) %*% lst[[k + 1]])
   G_pj[1, 2:ncol(G_pj)] <- as.vector(do.call(cbind, row_blocks))
   
   col_blocks <- lapply(1:(p + 1), function(k) lst[[k + 1]] %*% e)
   G_pj[2:nrow(G_pj), 1] <- as.vector(do.call(rbind, col_blocks))
   
   G_pj[2:nrow(G_pj), 2:ncol(G_pj)] <- toeplitz_block(p, lst)
   
   return(G_pj)
}

estimate_Gamma_h_list <- function(X, h_values, m = n_disc) {
   d <- ncol(X)
   lambda <- seq(0, 1, length.out = m)
   Gamma_h_list <- vector("list", length = length(h_values))
   n <- nrow(X)
   
   for (h in h_values) {
      mat <- array(NA_real_, dim = c(d, d, m))  # initialize with NA
      for (l in 1:m) {
         upper_index <- floor(lambda[l] * (n - h))
         if (upper_index > h) {  
            X1 <- X[1:upper_index, , drop = FALSE]
            X2 <- X[(1:upper_index) + h, , drop = FALSE]
            mat[, , l] <- t(X1) %*% (X2) / n
         }
      }
      Gamma_h_list[[h + 1]] <- mat
   }
   return(Gamma_h_list)
}

simulate_rejection_prob <- function(p, model = c("VAR", "LIN")) {
   model <- match.arg(model)
   h_values <- 0:p
   Gamma_list_true <- if (model == "VAR") Gamma_list_VAR else Gamma_list_LIN
   Sp <- compute_Sp_multi(p, Gamma_list_true)
   result_list <- list()
   
   for (N in sample_sizes) {
      rej_prob <- numeric(length(Delta))
      count_valid <- 0
      
      for (k in 1:rep) {
         eps_burn <- MASS::mvrnorm(100, rep(0, d), Sigma)
         eps <- MASS::mvrnorm(N + trunk, rep(0, d), Sigma)
         X_burn <- matrix(0, 100, d)
         X <- matrix(0, N, d)
         
         if (model == "VAR") {
            for (t in (q + 1):100) {
               X_burn[t, ] <- AR1 %*% X_burn[t - 1, ] + 
                  AR2 %*% X_burn[t - 2, ] + 
                  AR3 %*% X_burn[t - 3, ] + eps_burn[t, ]
            }
            X[1:q, ] <- X_burn[(100 - q + 1):100, ]
            for (t in (q + 1):N) {
               X[t, ] <- AR1 %*% X[t - 1, ] + 
                  AR2 %*% X[t - 2, ] + 
                  AR3 %*% X[t - 3, ] + eps[t, ]
            }
         } else {
            for (t in 1:N) {
               X[t, ] <- Reduce(`+`, lapply(0:trunk, function(j) {
                  Theta_list[[j + 1]] %*% eps[trunk + t - j, ]
               }))
            }
         }
         Gamma_array_list <- estimate_Gamma_h_list(X, h_values)
         S_p_hat_lambda <- rep(NA_real_, n_disc)
         for (l in 1:n_disc) {
            Gamma_list_lambda <- lapply(Gamma_array_list, function(arr) arr[, , l])
            val <- tryCatch(
               suppressWarnings(compute_Sp_multi(p, Gamma_list_lambda)),
               error = function(e) NA_real_
            )
            S_p_hat_lambda[l] <- val
         }
         if (is.na(S_p_hat_lambda[n_disc])) next
         valid_idx <- which(!is.na(S_p_hat_lambda))
         if (length(valid_idx) < 10) next  # require at least some usable lambdas
         
         V_N_S <- sum(lam_grid[valid_idx] * 
                         abs(S_p_hat_lambda[valid_idx] - S_p_hat_lambda[n_disc])) / n_disc
         Sp_hat <- S_p_hat_lambda[n_disc]
         
         rej_prob <- rej_prob + as.numeric(Sp_hat <= Delta + quantile(W_vals, alpha) * V_N_S)
         count_valid <- count_valid + 1
      }
      result_list[[as.character(N)]] <- if (count_valid > 0) rej_prob / count_valid else rep(NA_real_, length(Delta))
   }
   list(result = result_list, Sp = Sp)
}

rej_prob_results_VAR <- list()
rej_prob_results_LIN <- list()
Sp_VAR <- list()
Sp_LIN <- list()

out_VAR <- simulate_rejection_prob(p_value, model = "VAR")
out_LIN <- simulate_rejection_prob(p_value, model = "LIN")
  
rej_prob_results_VAR[[as.character(p_value)]] <- out_VAR$result
rej_prob_results_LIN[[as.character(p_value)]] <- out_LIN$result
Sp_VAR[[as.character(p_value)]] <- out_VAR$Sp
Sp_LIN[[as.character(p_value)]] <- out_LIN$Sp


##### Plots ####################################################################
par(mfrow = c(1, 2), oma = c(3, 1, 0.5, 1))

line_col   <- "black"
line_types <- c(3,4, 2, 1)
line_width <- 2

plot_measure_lines <- function(result_list, Sp, p, main_label, xlim_range, 
                               show_yaxis = TRUE
                      ) {
   plot(Delta, rep(0, length(Delta)), type = "n",
        xlim = xlim_range, ylim = c(0, 1),
        xlab = "", ylab = "", xaxt = "n", yaxt = "n",
        xaxs = "i", yaxs = "i"
   )
     
   title(main = main_label, line = 0.5, cex.main = 1.0, font.main = 1)
   
   abline(h = seq(0.1, 0.9, 0.1), col = "gray35", lty = "dotted")
   abline(v = seq(0, 1, 0.05), col = "gray35", lty = "dotted")
   abline(h = 0, col = "black")
   abline(h = alpha, col = "black", lty = 2)
   abline(h = 1, col = "black")
   if (!is.null(Sp) && !is.na(Sp)) {
      abline(v = Sp, col = "black", lty = 2)
      text(x = Sp, y = 0.80, labels = bquote(Delta == S[.(p)]),
      pos = 2, cex = 1.2)
   }
     
   for (i in seq_along(sample_sizes)) {
      y_vals <- result_list[[as.character(sample_sizes[i])]]
      if (!is.null(y_vals) && !all(is.na(y_vals))) {
         lines(Delta, y_vals, col = line_col, lwd = line_width, 
               lty = line_types[i]
         )
      }
   }
   axis(1, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1))
   if (show_yaxis) axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2))
}

par(mar = c(2, 2.5, 2, 2))
plot_measure_lines(rej_prob_results_VAR[[as.character(p_value)]],
                   Sp_VAR[[as.character(p_value)]],
                   p = p_value,
                   main_label = "VAR(3) Process",
                   xlim_range = c(0.2, 0.6),
                   show_yaxis = TRUE)

par(mar = c(2, 2, 2, 1))
plot_measure_lines(rej_prob_results_LIN[[as.character(p_value)]],
                   Sp_LIN[[as.character(p_value)]],
                   p = p_value,
                   main_label = "Multivariate Linear Process",
                   xlim_range = c(0.1, 0.8),
                   show_yaxis = FALSE)

par(xpd = NA)
legend(
  x = grconvertX(0.50, from = "ndc", to = "user"),
  y = grconvertY(0.13, from = "ndc", to = "user"),
  legend = paste("N =", sample_sizes),
  col = line_col, lty = line_types, lwd = line_width,
  horiz = TRUE, bty = "n", cex = 1.5, xjust = 0.5, x.intersp = 1,
  seg.len = 2.5, text.width = 0.175, adj = 0
)
par(xpd = FALSE)