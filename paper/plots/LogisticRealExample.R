#Display approximate and exact designs
#for the logistic regression model with two design variables

#Data in JZnoteLogistic.docx
#%% Reference paper: 
#% Haines, L. M., & Kabera, G. M. (2018). D-optimal designs 
#% for the two-variable binary logistic regression model 
#with interaction. Journal of Statistical Planning and Inference, 
#193, 136-150.
 

#criterion = "D";
#% beta = [1, 2, 2, 0.2]';
#beta = [-3, 4, 6, 1]' ; %Example 4.2 (b)
#S1 = [0, 1]; 
#S2 = [0, 1]; 
#p = 2; % Dimension



#design_app  
p1=c(         0,    0.2600,    0.1097)
p2=c(         0,    0.7400,    0.2470)
p3=c(    0.1600,    0.1400,    0.1416)
p4=c(    0.4000,         0,    0.0033)
p5=c(    0.6000,    0.4000,    0.2492)
p6=c(    1.0000,         0,    0.2492)

DesignApp=rbind(p1,p2,p3,p4,p5,p6)

x11=DesignApp[,1]
x21=DesignApp[,2]
w1=DesignApp[,3]

par(mfrow=c(2,2))
plot(x11,x21, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,1.1), ylim=c(0,1), col="red", pch=1,main="(a)")
text(x11+0.04,x21+0.05,w1,col="blue")

#From the paper (JSPI)
Dsupp1 <- c(1, 0, 0, 0.1758, 0.6049)
Dsupp2 <- c(0, 0.7351, 0.2649,0.1297, 0.4029)
ww1 <- c(0.2493, 0.2465, 0.1033, 0.1517, 0.2492)
points(Dsupp1, Dsupp2, col="green", pch = 6, cex = 2)
legend("topright", legend = c("Design 1", "Design 2"), 
       col = c("red","green"), pch = c(1, 6))



# exact design n= 10 ------------------------------------------------------
#design_ex =

g1=c(         0,    0.2695,    1)
g2=c(         0,    0.2707,    1)
g3=c(         0,    0.7460,    2)
g4=c(    0.4088,         0,    1)
g5=c(    0.5947,    0.3974,    1)
g6=c(    0.5966,    0.3971,    1)
g7=c(    0.5974,    0.3971,    1)
g8=c(    1.0000,         0,    2)

Dn10 <- rbind(g1,g2,g3,g4,g5,g6,g7,g8)
x12 <- Dn10[,1]
x22 <- Dn10[,2]
w2 <- Dn10[,3]
w2r <- c(2,2,1,3,2)

plot(x12, x22, xlab = expression(x[1]), ylab = expression(x[2]),
  xlim = c(0, 1.1), ylim = c(0, 1), col = "red", main = "(b)")
text(x12[c(1,3,4,6,8)] + 0.04, x22[c(1,3,4,6,8)] + 0.05, w2r, col="blue", cex = 1.1)


# n=15 --------------------------------------------------------------------

#design_ex =
h1=c(         0,    0.278,    1)
h2=c(         0,    0.2793,    1)
h3=c(         0,    0.7384,    1)
h4=c(         0,    0.7386,    1)
h5=c(         0,    0.7395,    1)
h6=c(     0.2300,    0.0907,    1)
h7=c(     0.2308,    0.0902,    1)
h8=c(     0.6001,    0.4020,    1)
h9=c(     0.6005,    0.4026,    1)
h10=c(    0.6009,    0.4035,    1)
h11=c(    0.6025,    0.4038,    1)
h12=c(    0.6026,    0.4028,    1)
h13=c(    1.0000,         0,    3)

Dn15=rbind(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13)
x13=Dn15[,1]
x23=Dn15[,2]
w3=Dn15[,3]
w3r=c(2,3, 2,5,3)

plot(x13, x23, xlab=expression(x[1]),ylab=expression(x[2]), 
     xlim=c(0,1.1), ylim=c(0,1), col="red", main="(c)")
text(x13[c(1,4,6,9,13)]+0.04,x23[c(1,4,6,9,13)]+0.05,w3r,col="blue")




# n=20 --------------------------------------------------------------------

#design_ex =
s1=c(         0,    0.263,    1)
s2=c(         0,    0.2632,    1)
s3=c(         0,    0.7348,    2)
s4=c(         0,    0.735,    1)
s5=c(         0,    0.7351,    1)
s6=c(         0,    0.7352,    1)
s7=c(         0.1684,    0.1342,    1)
s8=c(         0.1722,    0.1325,    1)
s9=c(    0.1737,    0.131,    1)
s10=c(    0.604,    0.4038,    1)
s11=c(    0.6047,    0.4024,    1)
s12=c(    0.6051,    0.4016,    1)
s13=c(    0.6062,    0.4035,    1)
s14=c(    0.6089,    0.402,    1)
s15=c(    1.0000,         0,    5)

Dn20=rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15)
x14=Dn20[,1]
x24=Dn20[,2]
w4=Dn20[,3]
w4r=c(2,5,3,5,5)

plot(x14,x24, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,1.1), ylim=c(0,1), col="red", main="(d)")
text(x14[c(1,3,7,10,15)]+0.04,x24[c(1,3,7,10,15)]+0.05,w4r,col="blue")



#Compute loss function value

LossD=function(Dmatrix, theta=c(-3,4,6,1))
{  p=dim(Dmatrix)[2]-1
   m=dim(Dmatrix)[1]
   wx=Dmatrix[,3]
   ww=wx/sum(wx)
   q=4
   Info=matrix(0,nrow=q,ncol=q)
   for (J in c(1:m))
   {  xx=Dmatrix[J,c(1:2)]
      fx=c(1,xx[1],xx[2],xx[1]*xx[2])
      u=sum(fx*theta)
      uw=exp(u)/(1+exp(u))^2
      Info=Info+ww[J]*uw*fx%*%t(fx)
    }
    1/(det(Info))^(1/q)
}

#Approximate designs

(L1=LossD(DesignApp)) #79.16625
Lapp/L1               #0.9998

#Exact designs
(L10=LossD(Dn10))    #80.48704
(L15=LossD(Dn15))    #80.90884

(L20=LossD(Dn20))    #79.16018

round(L1/L10,4)   #0.9836

round(L1/L15,4)   #0.9785
round(L1/L20,4)   # 1.0001

#D-opt in the paper
Dpaper=cbind(Dsupp1,Dsupp2,ww1)
(Lapp=LossD(Dpaper))  #79.15433


###################################
#Real Example
#################################
#D-optimal designs from the paper
#################################

#Xi-6
c1=c(0.0559,0)
c2=c(0.2689,0)
c3=c(0,0.3369)
c4=c(0,1.6193) 
c5=c(0.0261,0.1574) 
c6=c(0.1594,0.9603)
Xi6=rbind(c1,c2,c3,c4,c5,c6)
Xi6w=c(0.1145,0.2475,0.1145,0.2475,0.0260,0.25)
sum(Xi6w)
(Xi6=cbind(Xi6,Xi6w))

(D6L=LossD(Xi6, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #133.1108

#Xi-5
d1=c(0.0557,0) 
d2=c(0.2691,0) 
d3=c(0,0.3356) 
d4=c(0,1.6207) 
d5=c(0.1593,0.9596)
Xi5=rbind(d1,d2,d3,d4,d5)
Xi5w=c(0.1278, 0.2472,0.1278,0.2472,0.25)
sum(Xi5w)
(Xi5=cbind(Xi5,Xi5w))

(D5L=LossD(Xi5, theta=c(-2.2054,13.5803, 2.2547,1.6262)))
D6L/D5L

#Xi-4
e1=c(0.2665,0)
e2=c(0,1.6049)
e3=c(0.0260,0.1566)
e4=c(0.1603,0.9658)
Xi4=rbind(e1,e2,e3,e4)
Xi4w=c(0.25,0.25,0.25,0.25)
sum(Xi4w)
(Xi4=cbind(Xi4,Xi4w))

(D4L=LossD(Xi4, theta=c(-2.2054,13.5803, 2.2547,1.6262)))
D6L/D4L


#Optimal designs from CVX
#N=N1^2

#N1=21
a1=c(         0,    0.3000,    0.1502)
a2=c(         0,    1.6000,    0.2480)
a3=c(    0.1000,         0,    0.1150)
a4=c(    0.2000,    0.9000,    0.2500)
a5=c(    0.3000,         0,    0.2368)

DN21=rbind(a1,a2,a3,a4,a5)
(DN21L=LossD(DN21, theta=c(-2.2054,13.5803, 2.2547,1.6262)))  #136.9972

D6L/DN21L  #0.9716316

#Exact design, n=10
t1=c(         0,    0.3378,    1.0000)
t2=c(         0,    0.3508,    1.0000)
t3=c(         0,    1.6415,    1.0000)
t4=c(         0,    1.6466,    1.0000)
t5=c(    0.0733,         0,    1.0000)
t6=c(    0.1571,    0.9646,    1.0000)
t7=c(    0.1592,    0.9696,    1.0000)
t8=c(    0.1604,    0.9536,    1.0000)
t9=c(    0.2677,         0,    1.0000)
t10=c(    0.2690,         0,    1.0000)

DN21n10=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

(DN21n10L=LossD(DN21n10, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #135.5233
D6L/DN21n10L  #0.9821991

#N=31
#Approximate design

a1=c(         0,    0.3333,    0.1608)
a2=c(         0,    1.6000,    0.0883)
a3=c(         0,    1.6667,    0.1591)
a4=c(    0.0667,         0,    0.0953)
a5=c(    0.1333,    1.0000,    0.2344)
a6=c(    0.2000,    0.8667,    0.0162)
a7=c(    0.2667,         0,    0.2459)

DN31=rbind(a1,a2,a3,a4,a5,a6,a7)
(DN31L=LossD(DN31, theta=c(-2.2054,13.5803, 2.2547,1.6262)))  #134.4438 

D6L/DN31L  # 0.9900851

#Exact design

t1=c(         0,    0.3373,    1.0000)
t2=c(         0,    0.3469,    1.0000)
t3=c(         0,    1.6441,    1.0000)
t4=c(         0,    1.6461,    1.0000)
t5=c(    0.0741,         0,    1.0000)
t6=c(    0.1588,    0.9588,    1.0000)
t7=c(    0.1600,    0.9605,    1.0000)
t8=c(    0.1602,    0.9545,    1.0000)
t9=c(    0.2674,         0,    1.0000)
t10=c(    0.2691,         0,    1.0000)

DN31n10=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

(DN31n10L=LossD(DN31n10, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #135.52
 
D6L/DN31n10L  #0.9822227 


#N=41
#Approximate design

a1=c(         0,    0.3000,    0.2083)
a2=c(         0,    1.6500,    0.2489)
a3=c(    0.0500,         0,    0.0441)
a4=c(    0.1500,    1.0000,    0.2500)
a5=c(    0.2500,         0,    0.2488)

DN41=rbind(a1,a2,a3,a4,a5)
(DN41L=LossD(DN41, theta=c(-2.2054,13.5803, 2.2547,1.6262)))  #133.6368 

D6L/DN41L  # 0.9960644
  

#Exact design

t1=c(         0,    0.3384,    1.0000)
t2=c(         0,    0.3420,    1.0000)
t3=c(         0,    1.6410,    1.0000)
t4=c(         0,    1.6476,    1.0000)
t5=c(    0.0745,         0,    1.0000)
t6=c(    0.1585,    0.9597,    1.0000)
t7=c(    0.1590,    0.9621,    1.0000)
t8=c(    0.1599,    0.9615,    1.0000)
t9=c(    0.2681,         0,    1.0000)
t10=c(    0.2698,         0,    1.0000)

DN41n10=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

(DN41n10L=LossD(DN41n10, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #135.5199
 
 
D6L/DN41n10L  # 0.9822238


#N=51
#Approximate design

a1=c(         0,    0.3600,    0.0441)
a2=c(         0,    1.6000,    0.2487)
a3=c(    0.0400,    0.0800,    0.1985)
a4=c(    0.0800,         0,    0.0098)
a5=c(    0.1600,    0.9600,    0.2500)
a6=c(    0.2800,         0,    0.2489)

DN51=rbind(a1,a2,a3,a4,a5,a6)
(DN51L=LossD(DN51, theta=c(-2.2054,13.5803, 2.2547,1.6262)))  #133.3137
 

D6L/DN51L  #0.9984782

#Exact design


t1=c(         0,    0.3401,    1.0000)
t2=c(         0,    0.3435,    1.0000)
t3=c(         0,    1.6439,    1.0000)
t4=c(         0,    1.6443,    1.0000)
t5=c(    0.0745,         0,    1.0000)
t6=c(    0.1594,    0.9597,    1.0000)
t7=c(    0.1594,    0.9598,    1.0000)
t8=c(    0.1596,    0.9631,    1.0000)
t9=c(    0.2689,         0,    1.0000)
t10=c(    0.2690,         0,    1.0000)

DN51n10=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

(DN51n10L=LossD(DN51n10, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #135.5183 
 
 
D6L/DN51n10L  #0.982235  

#N=81
#Approximate design

a1=c(         0,    0.3600,    0.0441)
a2=c(         0,    1.6000,    0.2487)
a3=c(    0.0400,    0.0800,    0.1985)
a4=c(    0.0800,         0,    0.0098)
a5=c(    0.1600,    0.9600,    0.2500)
a6=c(    0.2800,         0,    0.2489)

DN81=rbind(a1,a2,a3,a4,a5,a6)
(DN81L=LossD(DN81, theta=c(-2.2054,13.5803, 2.2547,1.6262)))  #133.3137
 
 

D6L/DN81L  #0.9984782
 

#Exact design
t1=c(         0,    0.3401,    1.0000)
t2=c(         0,    0.3435,    1.0000)
t3=c(         0,    1.6439,    1.0000)
t4=c(         0,    1.6443,    1.0000)
t5=c(    0.0745,         0,    1.0000)
t6=c(    0.1594,    0.9597,    1.0000)
t7=c(    0.1594,    0.9598,    1.0000)
t8=c(    0.1596,    0.9631,    1.0000)
t9=c(    0.2689,         0,    1.0000)
t10=c(    0.2690,         0,    1.0000)
  
DN81n10=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)

(DN81n10L=LossD(DN81n10, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #135.5183
   
 
D6L/DN81n10L  #0.982235   


#N=51, n=15


t1=c(         0,    0.3669,    2.0000)
t2=c(         0,    1.6282,    1.0000)
t3=c(         0,    1.6287,    1.0000)
t4=c(         0,    1.6297,    1.0000)
t5=c(    0.0554,         0,    1.0000)
t6=c(    0.0555,         0,    1.0000)
t7=c(    0.1590,    0.9657,    1.0000)
t8=c(    0.1591,    0.9629,    1.0000)
t9=c(    0.1599,    0.9537,    1.0000)
t10=c(    0.1604,    0.9538,    1.0000)
t11=c(    0.2682,         0,    1.0000)
t12=c(    0.2691,         0,    2.0000)
t13=c(    0.2701,         0,    1.0000)


DN51n15=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13)

(DN51n15L=LossD(DN51n15, theta=c(-2.2054,13.5803, 2.2547,1.6262))) #133.9815
 
 
D6L/DN51n15L  #0.9935014

#Plot approximate and exact designs 

par(mfrow=c(2,2))

Xi5x1=Xi5[,1]
Xi5x2=Xi5[,2]

Xi6x1=Xi6[,1]
Xi6x2=Xi6[,2]

DN81x1=DN81[,1]
DN81x2=DN81[,2]
DN81w=DN81[,3]

DN51x1=DN51[,1]
DN51x2=DN51[,2]
DN51w=DN51[,3]

DN41x1=DN41[,1]
DN41x2=DN41[,2]
DN41w=DN41[,3]

par(mfrow=c(2,2))

plot(DN41x1,DN41x2, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,0.45), ylim=c(0,2), col="red", main="(a)")
points(Xi5x1,Xi5x2,col="green",pch=6,cex=1.5)
text(DN41x1+0.05,DN41x2+0.08,DN41w, col="blue")
legend("topright", legend=c("Design I", "Design II"), col=c("red","green"),
pch=c(1,6))


plot(DN51x1,DN51x2, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,0.45), ylim=c(0,2), col="red", main="(b)")
points(Xi6x1,Xi6x2,col="green",pch=6,cex=1.5)
text(DN51x1+0.05,DN51x2+0.08,DN51w, col="blue")
legend("topright", legend=c("Design III", "Design IV"), col=c("red","green"),
pch=c(1,6))



DN41n10x1=DN41n10[,1]
DN41n10x2=DN41n10[,2]
DN51n10x1=DN51n10[,1]
DN51n10x2=DN51n10[,2]
DN81n10x1=DN81n10[,1]
DN81n10x2=DN81n10[,2]

Dn10w=c(2,2,1,3,2)
#c(1,3,5,6,9)
plot(DN41n10x1,DN41n10x2, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,0.45), ylim=c(0,2), col="red", main="(c)")
text(DN41n10x1[c(1,3,5,6,9)]+0.03,DN41n10x2[c(1,3,5,6,9)]+0.08,Dn10w,col="blue")
points(DN81n10x1,DN81n10x2,col="green",pch=6,cex=1.5)
legend("topright", legend=c("Design V", "Design VI"), col=c("red","green"),
pch=c(1,6))

 
DN51n15x1=DN51n15[,1]
DN51n15x2=DN51n15[,2]
Dn15w=c(2,3,2,4,4)
#c(1,2,5,7,11)
plot(DN51n15x1,DN51n15x2, xlab=expression(x[1]),ylab=expression(x[2]),
xlim=c(0,0.45), ylim=c(0,2), col="red", main="(d)")
text(DN51n15x1[c(1,2,5,7,11)]+0.03,DN51n15x2[c(1,2,5,7,11)]+0.08,Dn15w,col="blue")