###################################################
# Calculation of the values of S_p in Section 3.4 #
###################################################

source("PivInf_AuxFunctions.R")

p_grid <- 1:7

# AR(5) process ################################################################
ar5_coefs <- c(-0.25, 0.10, 0.40, -0.25,  0.25) 
is_stationary(ar5_coefs) # Check for stationarity

# Calculate S_p values on p_grid ###############################################
gamma_AR5 <- true_ACV_uni_AR_via_ACFs(ar5_coefs, lag = max(p_grid))
S_vals <- sapply(p_grid, function(p) Sp_uni(gamma_AR5, p))
names(S_vals) <- p_grid
round(S_vals, 3)