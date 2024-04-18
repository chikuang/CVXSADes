#Approximate and exact maximin D-optimal designs

#%Compute maximin D-optimal design in Application 3
#%%Four dose-response models 
#%JCGS paper (Wong and Zhou)


#Nsim = 20;
#N = 201;         %number of design points   
#n = 10; %[10, 20, 30, 40, 50,60]';
#a = 0;  b = 500; %[a, b] is the design space
#u = linspace(a,b,N); %equally spaced N points in [a,b]
#v0 = 60; v1 = 294; v2=25;   %true parameter values for Emax I model 
#v02 = 60;  v12 = 340; v22 = 107.14; %true parameter values for Emax II model 
#v03 = 49.62; v13 = 290.51; v23 = 150; v33 = 45.51; %parameter values 

par(mfrow = c(2, 2), cex = 1)

#Approximate design

xa <- c(0, 20, 112.5, 205, 500)
wa <- c(0.241, 0.1789, 0.1314, 0.1248, 0.3239)

ya <- cumsum(wa)/sum(wa)

plot(xa, ya, xlab = "x", ylab = "distribution function", ylim = c(0, 1.2),
     xlim = c(-10, 515), col = "red", cex = 1.2, main = "(a)")
ma <- length(xa)


for (i in c(1:(ma - 1))) {
  lines(c(xa[i], xa[i + 1]), c(ya[i], ya[i]), col = "blue")
}

# Exact designs n=10

t1 <- c(0, 2)
t2 <- c(17.049, 1)
t3 <- c(17.062, 1)
t4 <- c(83.291, 1)
t5 <- c(111.48, 1)
t6 <- c(203.79, 1)
t7 <- c(500, 3)
Dn10 <- rbind(t1, t2, t3, t4, t5, t6, t7)
xn10 <- Dn10[, 1]
wn10 <- Dn10[, 2]
yn10 <- cumsum(wn10)/sum(wn10)

plot(xn10, yn10, xlab = "x", ylab = "distribution function", ylim = c(0,
                                                                      1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(b)")
m <- length(xn10)


for (i in c(1:(m - 1))) {
  lines(c(xn10[i], xn10[i + 1]), c(yn10[i], yn10[i]), col = "blue")
}

# n=20

t1 <- c(0, 6)
t2 <- c(18.035, 1)
t3 <- c(19.463, 1)
t4 <- c(33.542, 1)
t5 <- c(115.4, 1)
t6 <- c(115.49, 1)
t7 <- c(204.48, 1)
t8 <- c(204.63, 1)
t9 <- c(204.8, 1)
t10 <- c(500, 6)

Dn20 <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10)
xn20 <- Dn20[, 1]
wn20 <- Dn20[, 2]
yn20 <- cumsum(wn20)/sum(wn20)

plot(xn20, yn20, xlab = "x", ylab = "distribution function", 
     ylim = c(0, 1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(c)")
m <- length(xn20)


for (i in c(1:(m - 1))) {
  lines(c(xn20[i], xn20[i + 1]), c(yn20[i], yn20[i]), col = "blue")
}


# n=30
t2 <- c(14.2648,1)
t3 <- c(14.5472,1)
t4 <- c(14.9195,1)
t5 <- c(15.5030,1)
t6 <- c(18.4480,1)
t7 <- c(103.3939,1)
t8 <- c(105.1997,1)
t9 <- c(105.7545,1)
t10 <- c(106.2890,1)
t11 <- c(110.5053,1)
t12 <- c(202.1536,1)
t13 <- c(202.7992,1)
t14 <- c(203.5588,1)
t15 <- c(500.0000,9)

Dn30 <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14,
              t15)
xn30 <- Dn30[, 1]
wn30 <- Dn30[, 2]
yn30 <- cumsum(wn30)/sum(wn30)

plot(xn30, yn30, xlab = "x", ylab = "distribution function", ylim = c(0,
                                                                      1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(d)")
m <- length(xn30)


for (i in c(1:(m - 1))) {
  lines(c(xn30[i], xn30[i + 1]), c(yn30[i], yn30[i]), col = "blue")
}


for (i in c(1:(ma - 1))) {
  lines(c(xa[i], xa[i + 1]), c(ya[i], ya[i]), col = "green", lty = 2)
}

legend("bottomright", legend = c("approximate design", "exact design"),
       lty = c(2, 1), col = c("green", "blue"), cex = 0.75)
