#################################################
# Empirical rejection probabilities in Figure 1 #
#################################################

source("PivInf_AuxFunctions.R")
source("PivInf_GlobalParameters.R")

set.seed(123)

alpha    <- 0.05
p_values <- c(2, 4, 6)

# True values of S_p of linear processes #######################################
compute_S_vals <- function(filter_seq, a_name) {
   acv_vals <- true_ACV_uni_lin(filter_seq, max(p_values), sigma)
   S_vals <- sapply(p_values, function(p) Sp_uni(acv_vals, p))
   names(S_vals) <- paste0("S_", p_values)
   cat(a_name, ":", round(S_vals, 3), "\n")
   list(acv_vals = acv_vals, S_vals = S_vals)
}

a_poly <- 4
filter_seq_poly <- c(1,1,1,1,(4:trunk-2)^(-a_poly))
res_poly <- compute_S_vals(filter_seq_poly, "poly")

q_exp <- 0.85
filter_seq_exp <- c(2/3,2/3,2/3,q_exp^(3:trunk))
res_exp <- compute_S_vals(filter_seq_exp, "exp")

# Simulation depending on filter and ACVs ######################################
sim_rej_probs <- function(filter_seq, p) {
   rej_prob_list <- list()
   
   for (N in sample_sizes) {
      rej_count <- numeric(length(Delta))
      
      for (r in seq_len(rep)) {
         eps    <- rnorm(N + trunk, mean = 0, sd = sigma)
         X_full <- as.numeric(filter(eps, filter = filter_seq, sides = 1))
         X      <- X_full[(trunk + 1):(trunk + N)]
         
         h_values <- 0:p
         acv_h_result <- sapply(h_values, function(i) {
            sapply(lam_grid, function(lam) {
               upper_index <- floor(lam * (length(X) - i))
               if (upper_index > 0) {
                  sum(X[1:upper_index] * X[(1:upper_index) + i]) / length(X)
               } else 0
            })
         })
         
         S_p_hat_lam <- apply(acv_h_result, 1, Sp_uni, p = p)
         if (anyNA(S_p_hat_lam)) next
         
         V_N_S  <- mean(lam_grid * abs(S_p_hat_lam - S_p_hat_lam[n_disc]))
         Sp_hat <- Sp_uni(acv_h_result[n_disc, ], p)
         if (is.na(Sp_hat) || is.na(V_N_S)) next
         
         rej_count <- rej_count + as.numeric(
            Sp_hat <= Delta + quantile(W_vals, alpha) * V_N_S
         )
      }
      
      rej_prob_list[[as.character(N)]] <- rej_count / rep
   }
   rej_prob_list
}

rej_prob_results <- list(poly = list(), exp = list())
Sp_vals <- list(poly = list(), exp = list())

for (p in p_values) {
   rej_prob_results$poly[[as.character(p)]] <- sim_rej_probs(filter_seq_poly, p)
   rej_prob_results$exp [[as.character(p)]] <- sim_rej_probs(filter_seq_exp,  p)
   
   Sp_vals$poly[[as.character(p)]] <- Sp_uni(res_poly$acv_vals[1:(p + 1)], p)
   Sp_vals$exp [[as.character(p)]] <- Sp_uni(res_exp$acv_vals[1:(p + 1)],  p)
}

# Plot results #################################################################

layout(matrix(1:6, nrow = 3, byrow = TRUE))
par(oma = c(8, 4, 4, 3), mar = c(1.2, 2.5, 0.75, 0.5))

add_grid <- function(xlim_range, Sp) {
   abline(h = seq(0.1, 0.9, 0.1), col = "gray35", lty = "dotted")
   abline(v = seq(xlim_range[1], xlim_range[2], 0.05),
          col = "gray35", lty = "dotted")
   abline(h = c(0, 1), col = "black")
   abline(h = 0.05, col = "black", lty = 2, lwd = 1.4)
   abline(v = Sp, col = "black", lty = 2, lwd = 1.4)
}

plot_measure_lines <- function(result_list, Sp, p, row, col) {
   x_axis_on <- row == length(p_values)
   y_axis_on <- col == 1
   xlim_range <- if (col == 1) c(0.2, 0.8) else c(0.1, 0.6)
   
   plot(Delta, result_list[[as.character(sample_sizes[1])]], type = "n",
        xlim = xlim_range, ylim = c(0, 1),
        xlab = "", ylab = "", xaxt = "n", yaxt = "n",
        xaxs = "i", yaxs = "i"
   )
   
   add_grid(xlim_range, Sp)
   
   # NO panel letters (A–F) here
   
   text(Sp, 0.80, labels = bquote(Delta == S[.(p)]),
        pos = 2, cex = 1.3)
   
   for (i in seq_along(sample_sizes)) {
      lines(Delta, result_list[[as.character(sample_sizes[i])]],
            col = "black", lwd = 2,
            lty = c(3, 4, 2, 1)[i])
   }
   
   if (x_axis_on)
      axis(1, at = seq(xlim_range[1], xlim_range[2], 0.1))
   if (y_axis_on)
      axis(2, at = seq(0, 1, 0.1))
}

for (i in seq_along(p_values)) {
   p_str <- as.character(p_values[i])
   
   plot_measure_lines(
      rej_prob_results$poly[[p_str]],
      Sp_vals$poly[[p_str]],
      p_values[i], i, 1
   )
   
   plot_measure_lines(
      rej_prob_results$exp[[p_str]],
      Sp_vals$exp[[p_str]],
      p_values[i], i, 2
   )
}

mtext(bquote(theta[j] ~ 1/(j-2)^.(a_poly)),
      side = 3, outer = TRUE, at = 0.25, line = 1, cex = 1.1)
mtext(bquote(theta[j] ~ .(q_exp)^j),
      side = 3, outer = TRUE, at = 0.75, line = 1, cex = 1.1)

for (j in seq_along(p_values)) {
   mtext(bquote(p == .(p_values[j])),
         side = 2, outer = TRUE,
         at = (length(p_values) - j + 0.5) / length(p_values),
         line = 0.5, cex = 1)
}

par(xpd = NA)
legend(grconvertX(0.5, "ndc", "user"),
       grconvertY(0.055, "ndc", "user"),
       legend = paste("N =", sample_sizes),
       col = "black", lty = c(3, 4, 2, 1), lwd = 2,
       horiz = TRUE, bty = "n", cex = 2,
       xjust = 0.5, x.intersp = 1, seg.len = 2.5)
par(xpd = FALSE)