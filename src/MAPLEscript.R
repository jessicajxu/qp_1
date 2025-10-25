## load required libraries ========
library(here)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(magrittr)
library(data.table)
library(readr)
library(MAPLE)
library(tidyverse)

set.seed(268734) # set random seed

## read in summary data for exposure and outcome =====
exp = fread(here("data", "betax.assoc.txt"), head = T)
exp_raw = exp[,c("rs", "beta", "se", "af", "allele1", "allele0", "p_wald", "n_obs")] 
colnames(exp_raw) = c("SNP", "b", "se", "frq_A1", "A1", "A2", "P", "N")
out = fread(here("data","betay.assoc.txt"), head = T)
out_raw = out[,c("rs", "beta", "se", "af", "allele1", "allele0", "p_wald", "n_obs")] 
colnames(out_raw) = c("SNP", "b", "se", "frq_A1", "A1", "A2", "P", "N")

## estimate Omega (sample structure) using summary data and LD scores (LDSC) =====
paras = est_SS(dat1 = exp_raw,
               dat2 = out_raw,
               trait1.name = "exp",
               trait2.name = "out",
               ldscore.dir = here("data/eur_w_ld_chr"))

## read in z-scores for the candidate SNPs ======
zscorex = fread(here("data", "zscorex.txt"), head = F)
Zscore_1 = as.vector(zscorex[[1]])
zscorey = fread(here("data", "zscorey.txt"), head = F)
Zscore_2 = as.vector(zscorey[[1]])

## load the LD matrix =====
Sigmaxin = fread(here("data", "Sigmax.txt"), head = F)
Sigma1in = as.matrix(Sigmaxin)
Sigmayin = fread(here("data", "Sigmay.txt"), head = F)
Sigma2in = as.matrix(Sigmayin)

## sample sizes for exposure and outcome data
samplen1 = 20000
samplen2 = 20000

## estimates of Omega using LDSC
t1 = paras$Omega[1,1]
t2 = paras$Omega[2,2]
t12 = paras$Omega[1,2]

## run MAPLE ====
start1 <- Sys.time()
result = MAPLE(Zscore_1, Zscore_2, Sigma1in, Sigma2in, samplen1, samplen2,
               Gibbsnumber=1000, burninproportion=0.2, pi_beta_shape=0.5, pi_beta_scale=4.5,
               pi_c_shape=0.5, pi_c_scale=9.5, pi_1_shape=0.5, pi_1_scale=1.5,
               pi_0_shape=0.05, pi_0_scale=9.95, t1, t2, t12)
end1 <- Sys.time()
end1 - start1 # runtime

## Incorporating LDSC uncertainty through Omega simulations =========

t1_hat <- paras$Omega[1,1]
t2_hat <- paras$Omega[2,2]
t12_hat <- paras$Omega[1,2]
t1_se <- paras$Omega.se[1,1]
t2_se <- paras$Omega.se[2,2]
t12_se <- paras$Omega.se[1,2]

# simulating R samples of Omega
R <- 500
# Run MAPLE for each Omega sample
alpha_hat <- numeric(R)
alpha_se  <- numeric(R)
alpha_p <- numeric(R)

start2 <- Sys.time()
for (r in 1:R) {
  t1_r  <- rnorm(1, t1_hat, t1_se)
  t2_r  <- rnorm(1, t2_hat, t2_se)
  t12_r <- rnorm(1, t12_hat, t12_se)
  fit <- MAPLE(Zscore_1, Zscore_2, Sigma1in, Sigma2in, samplen1, samplen2,
              Gibbsnumber=1000, burninproportion=0.2, pi_beta_shape=0.5, pi_beta_scale=4.5,
              pi_c_shape=0.5, pi_c_scale=9.5, pi_1_shape=0.5, pi_1_scale=1.5,
              pi_0_shape=0.05, pi_0_scale=9.95, t1_r, t2_r, t12_r)
  alpha_hat[r] <- fit$causal_effect
  alpha_se[r]  <- fit$cause.se
  alpha_p[r] <- fit$causal_pvalue

  }
end2 <- Sys.time()
end2 - start2 # runtime

# results from original MAPLE
result$causal_effect
result$causal_pvalue
result$cause.se
# results from Omega simulated MAPLE
alpha_mean <- mean(alpha_hat)
alpha_total_var <- mean(alpha_se^2) + var(alpha_hat)
alpha_total_se <- sqrt(alpha_total_var)

# plotting
alpha_plot <- data.frame(alpha_hat = alpha_hat)

ggplot(alpha_plot, aes(x = alpha_hat)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black") +
  geom_vline(aes(xintercept = alpha_mean, color = "LDSC simulated"), linewidth = 1.2) +
  geom_vline(aes(xintercept = result$causal_effect, color = "Original"), linewidth = 1.2) +
  scale_color_manual(
    name = expression(hat(alpha)),
    values = c(`LDSC simulated` = "red", Original = "darkgreen")
  ) +
  labs(
    x = expression("Simulated"~hat(alpha)~"under LDSC uncertainty"),
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top"
  )
