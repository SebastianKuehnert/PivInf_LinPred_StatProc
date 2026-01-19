################################################################
# Coverage of PACs (lag=2,4) of AR(p) (p=2,4,6) in Section 4.3 #
################################################################

source("PivInf_GlobalParameters.R")
source("PivInf_AuxFunctions.R")

set.seed(123)

# Required parameters in Section 4.3 ###########################################
alpha <- 0.1                        # significance level
ar2_coefs <- c(-0.2, -0.3)          # is_stationary(ar2_coefs)
ar4_coefs <- c(ar2_coefs, 0.3, 0.2) # is_stationary(ar4_coefs)
ar6_coefs <- c(ar4_coefs, 0.1, 0.1) # is_stationary(ar6_coefs)

# Method BNS (execution of all may take ~ 5-6 min) #############################
emp_PACF_coverage_BNS(alpha, 2, ar2_coefs) # BNS for p=2 in scenario (i)
emp_PACF_coverage_BNS(alpha, 4, ar4_coefs) # BNS for p=4 in scenario (i)
emp_PACF_coverage_BNS(alpha, 2, ar6_coefs) # BNS for p=2 in scenario (ii)
emp_PACF_coverage_BNS(alpha, 4, ar6_coefs) # BNS for p=4 in scenario (ii)

# Method PIV (execution of all may take ~ 2-3 min) #############################
emp_PACF_coverage_PIV(alpha, 2, ar2_coefs) # PIV for p=2 in scenario (i)
emp_PACF_coverage_PIV(alpha, 4, ar4_coefs) # PIV for p=4 in scenario (i)
emp_PACF_coverage_PIV(alpha, 2, ar6_coefs) # PIV for p=2 in scenario (ii)
emp_PACF_coverage_PIV(alpha, 4, ar6_coefs) # PIV for p=4 in scenario (ii)