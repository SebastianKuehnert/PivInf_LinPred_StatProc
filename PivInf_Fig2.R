##########################
# Histograms in Figure 2 #
##########################

source("PivInf_AuxFunctions.R")
source("PivInf_GlobalParameters.R")

set.seed(123)

sample_sizes <- c(100, 200, 500, 1000, 5000)
alpha        <- 0.1
p_grid       <- 1:7
q_W          <- quantile(W_vals, alpha)
fill_cols    <- c("gray95", "gray75", "gray55", "gray35", "gray15")

# Plot function ################################################################
plot_hist_row <- function(p_hats_list,
                          show_titles = TRUE,
                          main_cex = 2,
                          title_line = 2) {
   
   for (j in seq_along(sample_sizes)) {
      N      <- sample_sizes[j]
      values <- na.omit(p_hats_list[[j]])
      values <- values[values >= 1 & values <= max(p_grid)]
      
      h <- hist(values,
                breaks = seq(0.5, max(p_grid) + 0.5, by = 1),
                plot   = FALSE)
      
      plot(0, 0, type = "n",
           xlim = range(h$breaks), ylim = c(0, 1),
           xlab = "", ylab = "", axes = FALSE)
      
      for (yy in seq(0, 1, 0.1)) {
         segments(x0 = 0.5, x1 = max(p_grid) + 0.5,
                  y0 = yy, y1 = yy,
                  col = "gray35", lty = "dotted")
      }
      
      for (i in seq_along(h$counts)) {
         rect(h$breaks[i], 0, h$breaks[i + 1], h$density[i],
              col = fill_cols[j], border = "black")
      }
      
      # column title (sample size) — only if requested (top row)
      if (show_titles) {
         title(main = bquote(N == .(N)), cex.main = main_cex, line = title_line)
      }
      
      axis(1, at = 1:max(p_grid), labels = 1:max(p_grid))
      if (j == 1) {
         axis(2, at = seq(0, 1, 0.1),
              labels = paste0(seq(0, 100, 10), "%"),
              las = 2)
      }
   }
   
   invisible(NULL)
}

# AR simulation ################################################################
nu_ar     <- 0.6
ar5_coefs <- c(-0.25, 0.10, 0.40, -0.25, 0.25)

compute_p_hats_ar <- function(N, p_grid, rep, ar_coefs, q_W, nu) {
   p_hats <- numeric(rep)
   for (r in 1:rep) {
      X <- arima.sim(model = list(ar = ar_coefs), n = N, sd = sigma)
      
      p_found <- NA
      for (p in p_grid) {
         Sp_seq <- sapply(1:n_disc, function(j) {
            acv_vec <- sapply(0:p, function(h) {
               upper <- floor(lam_grid[j] * (length(X) - h))
               if (upper > 0) {
                  sum(X[1:upper] * X[(1:upper) + h]) / length(X)
               } else {
                  NA
               }
            })
            if (anyNA(acv_vec)) return(NA)
            Sp_uni(matrix(acv_vec, ncol = 1), p)
         })
         
         Sp_seq <- na.omit(Sp_seq)
         if (length(Sp_seq) == 0) next
         
         Sp_hat <- Sp_seq[length(Sp_seq)]
         V_hat  <- sum(((1:length(Sp_seq)) / length(Sp_seq)) *
                          abs(Sp_seq - Sp_hat)) / n_disc
         
         if (!is.na(Sp_hat) && !is.na(V_hat) && Sp_hat < (1 - nu - q_W * V_hat)) {
            p_found <- p
            break
         }
      }
      p_hats[r] <- ifelse(is.na(p_found), NA, p_found)
   }
   p_hats
}

p_hats_ar_list <- lapply(sample_sizes, function(N)
   compute_p_hats_ar(N, p_grid, rep, ar5_coefs, q_W, nu_ar)
)

# Linear process simulation ####################################################
nu_lp      <- 0.65
filter_seq <- c(1, 1, 1, 1, (4:trunk - 2)^(-a))

compute_p_hats_lp <- function(N, p_grid, rep, q_W, nu, filter_seq) {
   p_hats <- numeric(rep)
   for (r in 1:rep) {
      
      eps    <- rnorm(N + trunk, mean = 0, sd = sigma)
      X_full <- as.numeric(stats::filter(eps, filter = filter_seq, sides = 1))
      X      <- X_full[(trunk + 1):(trunk + N)]
      
      p_found <- NA
      for (p in p_grid) {
         Sp_seq <- sapply(1:n_disc, function(j) {
            acv_vec <- sapply(0:p, function(h) {
               upper <- floor(lam_grid[j] * (length(X) - h))
               if (upper > 0) {
                  sum(X[1:upper] * X[(1:upper) + h]) / length(X)
               } else {
                  NA
               }
            })
            if (anyNA(acv_vec)) return(NA)
            Sp_uni(matrix(acv_vec, ncol = 1), p)
         })
         
         Sp_seq <- na.omit(Sp_seq)
         if (length(Sp_seq) == 0) next
         
         Sp_hat <- Sp_seq[length(Sp_seq)]
         V_hat  <- sum(((1:length(Sp_seq)) / length(Sp_seq)) *
                          abs(Sp_seq - Sp_hat)) / n_disc
         
         if (!is.na(Sp_hat) && !is.na(V_hat) && Sp_hat < (1 - nu - q_W * V_hat)) {
            p_found <- p
            break
         }
      }
      p_hats[r] <- ifelse(is.na(p_found), NA, p_found)
   }
   p_hats
}

p_hats_lp_list <- lapply(sample_sizes, function(N)
   compute_p_hats_lp(N, p_grid, rep, q_W, nu_lp, filter_seq)
)

# Plot #########################################################################
par(mfrow = c(2, length(sample_sizes)),
    oma  = c(1, 1, 0, 0), 
    mgp  = c(5, 0.5, 0),
    xpd  = NA)

par(mar = c(1, 1.7, 3.75, 0))
plot_hist_row(p_hats_ar_list,
              show_titles = TRUE,
              title_line = 2.75)

par(mar = c(1, 1.7, 3.75, 0))
plot_hist_row(p_hats_lp_list,
              show_titles = FALSE)