par(mfrow = c(2, 2), cex = 1)

xa <- c(0, 27.5, 30, 82.5, 85, 192.5, 500)
wa <- c(0.4621, 0.1155, 0.0084, 0.0735, 0.1168, 0.051, 0.1727)
sum(wa)

ya <- cumsum(wa)/sum(wa)


plot(xa, ya, xlab = "x", ylab = "distribution function", ylim = c(0, 1.2),
     xlim = c(-10, 515), col = "red", cex = 1.2, main = "(a)")
ma <- length(xa)


for (i in c(1:(ma - 1))) {
  lines(c(xa[i], xa[i + 1]), c(ya[i], ya[i]), col = "blue")
}


# Exact designs n=10 ------------------------------------------------------
t1 <- c(0, 4)
t2 <- c(11.2191, 1)
t3 <- c(37.2388, 1)
t4 <- c(89.1068, 1)
t5 <- c(183.974, 1)
t6 <- c(500, 2)



Dn10 <- rbind(t1, t2, t3, t4, t5, t6)
xn10 <- Dn10[, 1]
wn10 <- Dn10[, 2]
yn10 <- cumsum(wn10)/sum(wn10)

plot(xn10, yn10, xlab = "x", ylab = "distribution function", ylim = c(0,
                                                                      1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(b)")
m <- length(xn10)


for (i in c(1:(m - 1))) {
  lines(c(xn10[i], xn10[i + 1]), c(yn10[i], yn10[i]), col = "blue")
}



# n=20 --------------------------------------------------------------------
t1 <- c(0, 8)
t2 <- c(23.068, 1)
t3 <- c(28.471, 1)
t4 <- c(36.477, 1)

t5 <- c(45.771, 1)
t6 <- c(80.954, 1)
t7 <- c(83.54, 1)
t8 <- c(91.296, 1)
t9 <- c(173.32, 1)
t10 <- c(217.51, 1)
t11 <- c(500, 1)
t12 <- c(500, 3)

    
Dn20 <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11)
xn20 <- Dn20[, 1]
wn20 <- Dn20[, 2]
yn20 <- cumsum(wn20)/sum(wn20)

plot(xn20, yn20, xlab = "x", 
     ylab = "distribution function", 
     ylim = c(0, 1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(c)")
m <- length(xn20)


for (i in c(1:(m - 1))) {
  lines(c(xn20[i], xn20[i + 1]), c(yn20[i], yn20[i]), col = "blue")
}


# n=30 --------------------------------------------------------------------
t1 <- c(0, 13)
t2 <- c(7.6376, 1)
t3 <- c(26.7806, 1)
t4 <- c(27.5205, 1)
t5 <- c(38.7469, 1)
t6 <- c(54.6413, 1)
t7 <- c(74.8983, 1)
t8 <- c(75.7292, 1)
t9 <- c(83.3478, 1)
t10 <- c(87.5968, 1)
t11 <- c(89.3609, 1)
t12 <- c(189.4826, 1)
t13 <- c(191.4374, 1)
t14 <- c(500, 5)



Dn30 <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14)
xn30 <- Dn30[, 1]
wn30 <- Dn30[, 2]
yn30 <- cumsum(wn30)/sum(wn30)

plot(xn30, yn30, xlab = "x", ylab = "distribution function", 
     ylim = c(0, 1.2), xlim = c(-10, 515), col = "red", cex = 1.2, main = "(d)")
m <- length(xn30)


for (i in c(1:(m - 1))) {
  lines(c(xn30[i], xn30[i + 1]), c(yn30[i], yn30[i]), col = "blue")
}


for (i in c(1:(ma - 1))) {
  lines(c(xa[i], xa[i + 1]), c(ya[i], ya[i]), col = "green", lty = 2)
}

legend("bottomright", legend = c("approximate design", "exact design"),
       lty = c(2, 1), col = c("green", "blue"), cex = 0.75)