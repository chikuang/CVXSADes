# Display approximate and exact designs for the logistic regression
# model with two design variables

# Data in JZnoteLogistic.docx %% Reference paper: % Haines, L. M., &
# Kabera, G. M. (2018). D-optimal designs % for the two-variable
# binary logistic regression model with interaction. Journal of
# Statistical Planning and Inference, 193, 136-150.


# criterion = 'D'; % beta = [1, 2, 2, 0.2]'; beta = [-3, 4, 6, 1]' ;
# %Example 4.2 (b) S1 = [0, 1]; S2 = [0, 1]; p = 2; % Dimension


par(mfrow = c(2, 2), cex = 1)
# design_app

p1 <- c(0, 0.26, 0.1097)
p2 <- c(0, 0.74, 0.247)
p3 <- c(0.16, 0.14, 0.1416)
p4 <- c(0.4, 0, 0.0033)
p5 <- c(0.6, 0.4, 0.2492)
p6 <- c(1, 0, 0.2492)

DesignApp <- rbind(p1, p2, p3, p4, p5, p6)

x11 <- DesignApp[, 1]
x21 <- DesignApp[, 2]
w1 <- DesignApp[, 3]
plot(x11, x21, xlab = expression(x[1]), 
     ylab = expression(x[2]), xlim = c(0, 1.1), 
     ylim = c(0, 1), col = "red", pch = 1, main = "(a)")
text(x11 + 0.04, x21 + 0.05, w1, col = "blue")

# From the paper (JSPI)
Dsupp1 <- c(1, 0, 0, 0.1758, 0.6049)
Dsupp2 <- c(0, 0.7351, 0.2649, 0.1297, 0.4029)
ww1 <- c(0.2493, 0.2465, 0.1033, 0.1517, 0.2492)
points(Dsupp1, Dsupp2, col = "green", pch = 6, cex = 1.5)
legend("topright", legend = c("Design 1", "Design 2"), col = c("red", "green"),
       pch = c(1, 6))

# exact design n=10


# design_ex =

g1 <- c(0, 0.2706, 1)
g2 <- c(0, 0.2708, 1)
g3 <- c(0, 0.7458, 1)
g4 <- c(0, 0.7462, 1)
g5 <- c(0.4086, 0, 1)
g6 <- c(0.594, 0.3972, 1)
g7 <- c(0.5954, 0.3981, 1)
g8 <- c(0.598, 0.3978, 1)
g9 <- c(1, 0, 2)

Dn10 <- rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9)
x12 <- Dn10[, 1]
x22 <- Dn10[, 2]
w2 <- Dn10[, 3]
w2r <- c(2, 2, 1, 3, 2)

plot(x12, x22, xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0,
                                                                          1.1), ylim = c(0, 1), col = "red", main = "(b)")
text(x12[c(1, 3, 5, 6, 9)] + 0.04, x22[c(1, 3, 5, 6, 9)] + 0.05, w2r, col = "blue")

# exact design n=15


# design_ex =
h1 <- c(0, 0.2587, 1)
h2 <- c(0, 0.2597, 1)
h3 <- c(0, 0.7367, 1)
h4 <- c(0, 0.737, 1)
h5 <- c(0, 0.7373, 1)
h6 <- c(0, 0.7379, 1)
h7 <- c(0.195, 0.1188, 1)
h8 <- c(0.2008, 0.1135, 1)
h9 <- c(0.6007, 0.4025, 1)
h10 <- c(0.6039, 0.403, 1)
h11 <- c(0.6043, 0.4016, 1)
h12 <- c(0.6047, 0.4016, 1)
h13 <- c(1, 0, 3)

Dn15 <- rbind(h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13)
x13 <- Dn15[, 1]
x23 <- Dn15[, 2]
w3 <- Dn15[, 3]
w3r <- c(2, 4, 2, 4, 3)

plot(x13, x23, xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0,
                                                                          1.1), ylim = c(0, 1), col = "red", main = "(c)")
text(x13[c(1, 3, 7, 9, 13)] + 0.04, x23[c(1, 3, 7, 9, 13)] + 0.05, w3r,
     col = "blue")


# exact design n=20


# design_ex =
s1 <- c(0, 0.2573, 1)
s2 <- c(0, 0.2576, 1)
s3 <- c(0, 0.2583, 1)
s4 <- c(0, 0.7413, 1)
s5 <- c(0, 0.742, 1)
s6 <- c(0, 0.7422, 1)
s7 <- c(0, 0.7426, 1)
s8 <- c(0, 0.7427, 1)
s9 <- c(0.2007, 0.1142, 1)
s10 <- c(0.2034, 0.1141, 1)
s11 <- c(0.601, 0.4012, 1)
s12 <- c(0.6013, 0.4018, 1)
s13 <- c(0.6017, 0.4027, 1)
s14 <- c(0.6031, 0.4007, 1)
s15 <- c(0.6038, 0.4006, 1)
s16 <- c(1, 0, 5)

Dn20 <- rbind(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14,
              s15, s16)
x14 <- Dn20[, 1]
x24 <- Dn20[, 2]
w4 <- Dn20[, 3]
w4r <- c(3, 5, 2, 5, 5)

plot(x14, x24, xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0,
                                                                          1.1), ylim = c(0, 1), col = "red", main = "(d)")
text(x14[c(1, 4, 9, 11, 16)] + 0.04, x24[c(1, 4, 9, 11, 16)] + 0.05, w4r,
     col = "blue")



# Compute loss function value

LossD <- function(Dmatrix, theta = c(-3, 4, 6, 1)) {
  p <- dim(Dmatrix)[2] - 1
  m <- dim(Dmatrix)[1]
  wx <- Dmatrix[, 3]
  ww <- wx/sum(wx)
  q <- 4
  Info <- matrix(0, nrow = q, ncol = q)
  for (J in c(1:m)) {
    xx <- Dmatrix[J, c(1:2)]
    fx <- c(1, xx[1], xx[2], xx[1] * xx[2])
    u <- sum(fx * theta)
    uw <- exp(u)/(1 + exp(u))^2
    Info <- Info + ww[J] * uw * fx %*% t(fx)
  }
  1/(det(Info))^(1/q)
}

# Approximate designs

(L1 <- LossD(DesignApp))  #79.16625
# Lapp/L1  #0.9998

# Exact designs
(L10 <- LossD(Dn10))  #80.48721
(L15 <- LossD(Dn15))  #79.73108

(L20 <- LossD(Dn20))  #79.18744

round(L1/L10, 4)  #0.9836

round(L1/L15, 4)  #0.9929
round(L1/L20, 4)  # 0.9997

# D-opt in the paper
Dpaper <- cbind(Dsupp1, Dsupp2, ww1)
(Lapp <- LossD(Dpaper))  #79.15433
