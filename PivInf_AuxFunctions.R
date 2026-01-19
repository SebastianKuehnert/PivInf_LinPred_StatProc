##################################################
# Auxiliary functions throughout all simulations #
##################################################

library(matrixcalc)
source("PivInf_GlobalParameters.R")

########## I: General functions ################################################

# Testing stationarity of AR process given vector of AR coefficients
is_stationary <- function(ar_coefs) {
   roots <- polyroot(c(1, - ar_coefs))
   return(all(Mod(roots) > 1))
}

# Quantiles of W in Eq. (2.12)
n_sim <- 10000
W_vals <- replicate(n_sim, {
   B <- cumsum(rnorm(n_disc)) / sqrt(n_disc)
   B1 <- B[n_disc]
   B1 * n_disc / sum(abs(B[1:n_disc] - (1:n_disc / n_disc) * B1))
})


########## II: Specific functions in the univariate setting ####################

# Univariate M_p in Eq. (1.1)
Mp_uni <- function(gamma, p) {
   if (!is.numeric(gamma)) stop("gamma must be numeric!")
   if (p < 0 || p >= length(gamma)) stop("p must be in [0, length(gamma)-1]!")
   if (p == 0) return(gamma[1])
   num <- det(toeplitz(gamma[1:(p + 1)]))
   denom <- det(toeplitz(gamma[1:p]))
   if (is.na(denom) || denom == 0) return(NA_real_)
   return(num / denom)
}

# Univariate S_p in Eq. (1.2)
Sp_uni <- function(gamma, p) {
   if (!is.numeric(gamma)) stop("gamma must be numeric!")
   if (p < 0 || p >= length(gamma)) stop("p must be in [0, length(gamma)-1]!")
   if (p == 0) return(1)
   num <- Mp_uni(gamma, p)
   denom <- Mp_uni(gamma, 0)
   if (is.na(denom) || denom == 0) return(NA_real_)
   return(num / denom)
}

# Simulations of n realizations of AR process
simulate_ar <- function(ar_coefs, n, burn = 200, innov_sd = sigma) {
   x <- numeric(n + burn)
   eps <- rnorm(n + burn, 0, innov_sd)
   for (t in (length(ar_coefs) + 1):(n + burn)) {
      x[t] <- sum(ar_coefs * x[(t - 1):(t - length(ar_coefs))]) + eps[t]
   }
   return(x[(burn + 1):(burn + n)])
}

# True ACVs of univariate linear processes #################################
true_ACV_uni_lin <- function(x, h_max, innov_sd = sigma){
   if (!is.numeric(x)) stop("x must be numeric!")
   if (h_max <0 || h_max > length(x)) stop("h_max has to be in [0,length(x)]!")
   gamma_vals <- numeric()
   for (h in 0:h_max) {
      upper <- length(x) - h
      gamma_vals[h + 1] <- innov_sd^2 * sum(x[1:upper] * x[(1:upper) + h])
   }
   return(gamma_vals)
}

# True ACVs of univariate AR processes via ACFs
true_ACV_uni_AR_via_ACFs <- function(ar_coefs, lag = length(ar_coefs), 
                                     innov_sd = sigma) {
   rho_vals <- ARMAacf(ar = ar_coefs, lag.max = lag)
   gamma0 <- innov_sd^2 / (1 - sum(ar_coefs * rho_vals[2:(length(ar_coefs)+1)]))
   return(gamma0 * rho_vals)
}

# True PACFs of univariate AR processes via ACFs

# Calculation of the PACF kappa_index of AR process with coefficients ar_coefs
true_PACF_uni_AR <- function(index, ar_coefs, innov_sd = sigma) {
   acv <- true_ACV_uni_AR_via_ACFs(ar_coefs, length(ar_coefs), innov_sd)
   G <- toeplitz(acv[1:index])
   rhs <- acv[2:(index + 1)]
   if (!is.positive.definite(G)) return(NA_real_)
   return(tail(solve(G, rhs), 1))
}

# Empirical sequential ACVs of univariate processes
emp_ACV <- function(x, h, lam) {
   n <- length(x) # sample size
   m <- floor(lam * (n - h))
   if (m <= 0) return(NA_real_)
   return(sum(x[1:m]*x[(1:m) + h]) / n)
}

# Empirical sequential rth PACF of univariate AR process
emp_PACF_uni_AR <- function(index, x, lam) {
   gammas <- sapply(0:index, function(h) emp_ACV(x, h, lam))
   if (any(is.na(gammas))) return(NA_real_)
   G   <- toeplitz(gammas[1:index])
   rhs <- gammas[2:(index + 1)]
   if (!is.positive.definite(G)) return(NA_real_)
   return(tail(solve(G, rhs), 1))
}

# Empirical coverage via method BNS
# The following function calculates the empirical coverage for the PACF 
# \kappa_index to significance level alpha with the BNS method for a given AR 
# process determined by its coefficients ar_coefs (for several sample sizes)
emp_PACF_coverage_BNS <- function(alpha, index, ar_coefs, innov_sd = sigma) {
   true_kappa <- true_PACF_uni_AR(index, ar_coefs)
   
   results <- data.frame()
   
   for (N in sample_sizes) {
      coverage <- rep(NA_real_, rep)
      length   <- rep(NA_real_, rep)
      
      for (r in 1:rep) {
         x <- simulate_ar(ar_coefs, N, innov_sd = innov_sd)
         
         # fit AR(p) model
         fit <- tryCatch(
            arima(x, order = c(index, 0, 0), method = "ML"),
            error = function(e) NULL
         )
         if (is.null(fit)) next
         
         coef_name <- paste0("ar", index)
         est <- tryCatch(coef(fit)[coef_name], error = function(e) NA)
         se  <- tryCatch(sqrt(diag(vcov(fit)))[coef_name], error = function(e) NA)
         
         if (is.na(est) || is.na(se)) next
         
         ci_lower <- est - qnorm(1 - alpha/2) * se
         ci_upper <- est + qnorm(1 - alpha/2) * se
         
         coverage[r] <- (ci_lower <= true_kappa && ci_upper >= true_kappa)
         length[r]   <- ci_upper - ci_lower
      }
      
      results <- rbind(results, data.frame(
         N         = N,
         Coverage  = round(mean(coverage, na.rm = TRUE), 3),
         Length    = mean(length, na.rm = TRUE),
         ValidReps = sum(!is.na(coverage))
      ))
   }
   
   cat("\nEmp. rej. probs with BNS for PACF kappa_", index, 
       " of AR(", length(ar_coefs), ") process:\n", sep = "")
   print(results, row.names = FALSE)
   
   invisible(results)
}

##### Empirical coverage via method PIV ########################################
# The following function calculates the empirical coverage for the PACF 
# \kappa_index to significance level alpha with the method PIV for a given AR 
# process determined by its coefficients ar_coefs (for several sample sizes)
emp_PACF_coverage_PIV <- function(alpha, index, ar_coefs, innov_sd = sigma) {
   quant <- quantile(W_vals, 1 - alpha/2)
   true_kappa <- true_PACF_uni_AR(index, ar_coefs, innov_sd)
   
   results <- data.frame()
   
   for (N in sample_sizes) {
      coverage <- rep(NA_real_, rep)
      length   <- rep(NA_real_, rep)
      
      for (r in 1:rep) {
         x <- simulate_ar(ar_coefs, N, innov_sd = innov_sd)
         
         kappas_hat_seq <- sapply(lam_grid, function(lam) 
            emp_PACF_uni_AR(index, x, lam)
         )
         if (any(is.na(kappas_hat_seq))) next
         kappa_hat <- kappas_hat_seq[n_disc]
         
         V_hat <- sum(lam_grid * abs(kappas_hat_seq - kappa_hat)) / n_disc
         
         ci_lower <- kappa_hat - quant * V_hat
         ci_upper <- kappa_hat + quant * V_hat
         
         coverage[r] <- (ci_lower <= true_kappa && ci_upper >= true_kappa)
         length[r]   <- ci_upper - ci_lower
      }
      
      results <- rbind(results, data.frame(
         N         = N,
         Coverage  = round(mean(coverage, na.rm = TRUE), 3),
         Length    = mean(length,   na.rm = TRUE),
         ValidReps = sum(!is.na(coverage))
      ))
   }
   
   cat("\nEmp. rej. probs with PIV for PACF kappa_", index, 
       " of AR(", length(ar_coefs), ") process:\n", sep = "")
   print(results, row.names = FALSE)
   
   invisible(results)
}


########## III: Specific functions in the multivariate setting #################

##### Function that creates block toeplitz matrix ##############################
toeplitz_block <- function(p, lst) {
   if (p >= length(lst)) stop("There are not enough elements in the list to 
                              create the block toeplitz matrix!")
   if (p < 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
   if (p == 0) return(lst[[1]])
   d <- nrow(lst[[1]])
   T_p <- matrix(0, nrow = d * (p + 1), ncol = d * (p + 1))
   for (i in 0:p) {
      for (j in 0:p) {
         block <- if (j - i >= 0) {
            lst[[j - i + 1]]
         } else {
         t(lst[[abs(j - i) + 1]])
         }
         T_p[(i*d + 1):((i + 1)*d), (j*d + 1):((j + 1)*d)] <- block
      }
   }
   return(T_p)
}

##### M_p in Section 5 #########################################################
compute_Mp_multi <- function(p, lst) {
   d <- nrow(lst[[1]])
   
   if (p == 0) {
      return(sum(diag(lst[[1]])))
   } else {
      det_Gpmin1 <- det(toeplitz_block(p - 1, lst))
      if (!is.finite(det_Gpmin1) || abs(det_Gpmin1) < 1e-9) return(NA)  
      
      det_sum <- 0
      for (j in 1:d) {
         det_Gpj <- det(build_Gpj(p - 1, j, lst))
         if (!is.finite(det_Gpj)) return(NA)
         det_sum <- det_sum + det_Gpj
      }
      
      M_p <- det_sum / det_Gpmin1
      return(M_p)
   }
}

##### S_p in Section 5 #########################################################
compute_Sp_multi <- function(p, lst) {
   M_p <- compute_Mp_multi(p, lst)
   M_0 <- compute_Mp_multi(0, lst)  # = tr(Gamma_0)
   return(M_p / M_0)
}