# Application 2 from Haines et al. (2018)


design_opt <- t(matrix(c(0.0559, 0.2689, 0, 0, 0.0261, 0.1594,
                       0, 0, 0.3369, 1.6193, 0.1574, 0.9603,
                       0.1145, 0.2475, 0.1145, 0.2475, 0.0260, 0.25), nrow = 3,
                     byrow = T))

beta <- c(-2.2054, 13.5803, 2.2547, 1.6262)


calc_info_ex <- function(design_ex, beta){
  n <- sum(design_ex[,3])
  w <- design_ex[,3]/n
  n_pt <- nrow(design_ex)
  M <- matrix(0, nrow = 4, ncol = 4)
  for(i in 1:n_pt){
    xx <- design_ex[i, 1:2]
    rx <- c(1, xx, xx[1] * xx[2])
    Gamma <- as.numeric(exp(beta %*% rx)/(1+exp(beta %*% rx))^2)
    M <- M + w[i] * Gamma * tcrossprod(rx)
  }
  return(M)
}

calc_info_app <- function(design_app, beta){
  w <- design_app[,3]
  n_pt <- nrow(design_app)
  M <- matrix(0, nrow = 4, ncol = 4)
  for(i in 1:n_pt){
    xx <- design_app[i, 1:2]
    rx <- c(1, xx, xx[1] * xx[2])
    Gamma <- as.numeric(exp(beta %*% rx)/(1+exp(beta %*% rx))^2)
    M <- M + w[i] * Gamma * tcrossprod(rx)
  }
  return(M)
}

calc_D_opt <- function(M, q = length(beta)){
  exp(-log(det(M)^(1/q)))
}

# N =21 -------------------------------------------------------------------
design_app_N21 <- t(matrix(c(0, 0, 0.1, 0.2, 0.3,
                           0.3, 1.6, 0, 0.9, 0,
                           0.15017, 0.24804, 0.11501, 0.25, 0.23679), nrow = 3,
                    byrow = T))
design_ex_N21 <- matrix(c(0, 0.3436, 1,
                      0, 0.345, 1,
                      0, 1.6457, 1,
                      0, 1.6466, 1,
                      0.074, 0, 1,
                      0.1592, 0.9568, 1,
                      0.1596, 0.9575, 1,
                      0.1606, 0.9568, 1,
                      0.2681, 0, 1,
                      0.2685, 0, 1),
                    ncol = 3, byrow = TRUE)
Mopt <- calc_info_app(design_opt, beta)
loss_opt <- calc_D_opt(Mopt)
M1 <- calc_info_app(design_app_N21, beta)
loss_app_N21 <- calc_D_opt(M1)
M2 <- calc_info_ex(design_ex_N21, beta)
loss_ex_N21 <- calc_D_opt(M2)
eff_N_N21 <- loss_opt/loss_app_N21
eff_n_N21 <- loss_opt/loss_ex_N21
# N31 ---------------------------------------------------------------------

design_app_N31 <- matrix(c(0, 0.33333, 0.16084,
                       0, 1.6, 0.088287,
                       0, 1.6667, 0.15908,
                       0.066667, 0, 0.09531,
                       0.13333, 1, 0.23436,
                       0.2, 0.86667, 0.016209,
                       0.26667, 0, 0.24592),
                     ncol = 3, byrow = TRUE)
design_ex_N31 <- matrix(c(0, 0.3333, 1,
                      0, 1.6183, 1,
                      0, 1.6206, 1,
                      0.0556, 0, 1,
                      0.1579, 0.9623, 1,
                      0.1585, 0.9586, 1,
                      0.1592, 0.9597, 1,
                      0.1599, 0.9599, 1,
                      0.269, 0, 1,
                      0.2691, 0, 1),
                    ncol = 3, byrow = TRUE)
Mopt <- calc_info_app(design_opt, beta)
loss_opt <- calc_D_opt(Mopt)
M1_N31 <- calc_info_app(design_app_N31, beta)
loss_app_N31 <- calc_D_opt(M1_N31)
M2_N31 <- calc_info_ex(design_ex_N31, beta)
loss_ex_N31 <- calc_D_opt(M2_N31)
eff_N_N31 <- loss_opt/loss_app_N31
eff_n_N31 <- loss_opt/loss_ex_N31

# N41 ---------------------------------------------------------------------
design_app_N41 <- matrix(c(0, 0.3, 0.20828,
                       0, 1.65, 0.24886,
                       0.05, 0, 0.044061,
                       0.15, 1, 0.25,
                       0.25, 0, 0.2488),
                     ncol = 3, byrow = TRUE)

design_ex_N41 <- matrix(c(0, 1.5894, 1,
                      0, 1.6056, 1,
                      0.0241, 0.174, 1,
                      0.0328, 0.1064, 1,
                      0.161, 0.971, 1,
                      0.1611, 0.9649, 1,
                      0.1613, 0.9558, 1,
                      0.2662, 0, 1,
                      0.2666, 0, 1,
                      0.2695, 0, 1),
                    ncol = 3, byrow = TRUE)
Mopt <- calc_info_app(design_opt, beta)
loss_opt <- calc_D_opt(Mopt)
M1_N41 <- calc_info_app(design_app_N41, beta)
loss_app_N41 <- calc_D_opt(M1_N41)
M2_N41 <- calc_info_ex(design_ex_N41, beta)
loss_ex_N41 <- calc_D_opt(M2_N41)
eff_N_N41 <- loss_opt/loss_app_N41
eff_n_N41 <- loss_opt/loss_ex_N41

# N51 ---------------------------------------------------------------------
design_app_N51 <- matrix(c(0, 0.36, 0.044132,
                       0, 1.6, 0.24865,
                       0.04, 0.08, 0.19853,
                       0.08, 0, 0.0098219,
                       0.16, 0.96, 0.24997,
                       0.28, 0, 0.2489),
                     ncol = 3, byrow = TRUE)

design_ex_N51 <- matrix(c(0, 1.6037, 1,
                      0, 1.6048, 1,
                      0.0229, 0.1685, 1,
                      0.0291, 0.1362, 1,
                      0.1594, 0.9659, 1,
                      0.16, 0.9666, 1,
                      0.1613, 0.9669, 1,
                      0.2658, 0, 1,
                      0.2681, 0, 1,
                      0.269, 0, 1),
                    ncol = 3, byrow = TRUE)
Mopt <- calc_info_app(design_opt, beta)
loss_opt <- calc_D_opt(Mopt)
M1_N51 <- calc_info_app(design_app_N51, beta)
loss_app_N51 <- calc_D_opt(M1_N51)
M2_N51 <- calc_info_ex(design_ex_N51, beta)
loss_ex_N51 <- calc_D_opt(M2_N51)
eff_N_N51 <- loss_opt/loss_app_N51
eff_n_N51 <- loss_opt/loss_ex_N51

# N81 ---------------------------------------------------------------------

design_app_N81 <- matrix(c(0, 0.35, 0.086135,
                       0, 1.6, 0.24754,
                       0.05, 0.025, 0.16814,
                       0.15, 0.975, 0.24999,
                       0.275, 0, 0.24819),
                     ncol = 3, byrow = TRUE)

design_ex_N81 <- matrix(c(0, 0.4464, 1,
                      0, 1.6173, 1,
                      0, 1.6188, 1,
                      0.0569, 0, 1,
                      0.0573, 0, 1,
                      0.1591, 0.9599, 1,
                      0.1593, 0.9587, 1,
                      0.1595, 0.9589, 1,
                      0.2732, 0, 1,
                      0.2734, 0, 1),
                    ncol = 3, byrow = TRUE)
Mopt <- calc_info_app(design_opt, beta)
loss_opt <- calc_D_opt(Mopt)
M1_N81 <- calc_info_app(design_app_N81, beta)
loss_app_N81 <- calc_D_opt(M1_N81)
M2_N81 <- calc_info_ex(design_ex_N81, beta)
loss_ex_N81 <- calc_D_opt(M2_N81)
eff_N_N81 <- loss_opt/loss_app_N81
eff_n_N81 <- loss_opt/loss_ex_N81

eff <- matrix(c(eff_N_N21, eff_N_N31, eff_N_N41, eff_N_N51, eff_N_N81,
         eff_n_N21, eff_n_N31, eff_n_N41, eff_n_N51, eff_n_N81,
         19, 94, 12, 87, 51), nrow = 3, byrow =T)
rownames(eff) <- c("Eff_app", "Eff_ex", "seed")
round(eff, 4)
