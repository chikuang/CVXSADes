library(tibble)
library(dplyr)
partial_deri <- function(x, w){
  w * (c(1,x,x^2,x^3) %o% c(1,x,x^2,x^3))
}
calc_loss <- function(design, p){
  point <- design$x
  weight <- design$weight
  n_x <- nrow(design) 
  FIM <- matrix(0, nrow = p+1, ncol = p+1)
  for(i in 1:n_x){
    # print(partial_deri(point[i], weight[i]))
    FIM <- FIM + partial_deri(point[i], weight[i])
  }
  # -log(det(FIM)^(1/(p+1)))
  -log(det(FIM))
}

# N21 n20 -----------------------------------------------------------------
N <- 21; n <- 20; p <- 3; # cubic, so p = 3,  with #seed 244
design_app <- matrix(c(-1, -0.5, -0.4, 0.4  ,0.5, 1, 
                       0.2495,    0.1129,    0.1376,    0.1376,    0.1129,    0.2495),
                     byrow = T, nrow = 2)
design_ex <- matrix(c(-1, -0.4476, -0.4467, -0.4465, -0.445, 0.4452, 0.4467, 0.447, 0.4479, 0.4481, 1, 
                      5, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5),
                    byrow=T, nrow = 2)

design_weight <- design_ex[2,]/n
res_design_app <- tibble(x = design_app[1,], n_point = NA, weight = design_app[2,])
res_design_ex <- tibble(x = design_ex[1,], n_point = design_ex[2,], weight = design_weight)

calc_loss(res_design_app, 3) #5.28968
calc_loss(res_design_ex, 3) #5.274609
# loss_appr <- 1.3224
# loss_ex <- 1.3187


# N21m n18 ----------------------------------------------------------------
N <- 21; n <- 18;  #seed 186
design_app <- matrix(c(-1, -0.5, -0.4, 0.4  ,0.5, 1, 
                       0.2495,    0.1129,    0.1376,    0.1376,    0.1129,    0.2495),
                     byrow = T, nrow = 2)


design_ex <- matrix(c(-1, -0.4477, -0.4469,  -0.4467, -0.4464, -0.4458,  0.446, 0.4475, 0.4476, 0.4483, 1, 
                      4, 1, 1, 1, 1, 1, 1, 2, 1,1,4),
                    byrow = T, nrow = 2)
design_weight <- design_ex[2,]/n
res_design_app <- tibble(x = design_app[1,], n_point = NA, weight = design_app[2,])
res_design_ex <- tibble(x = design_ex[1,], n_point = design_ex[2,], weight = design_weight)

calc_loss(res_design_app, 3) #5.28968
calc_loss(res_design_ex, 3)  #5.29945
# loss_appr <- 1.3224
# loss_ex <- 1.3249

# N 51, n 20 --------------------------------------------------------------
N <- 51; n <- 20; # seed 435
design_app <- matrix(c(-1, -0.44, 0.44, 1,
                       0.25, 0.25, 0.25, 0.25),
                     byrow=T,nrow=2)

design_ex <- matrix(c(-1, -0.4482,  -0.4476, -0.4474, -0.447,  -0.4459, 0.4464, 0.4469, 0.4474, 0.4476, 1,
                      5, 1, 1, 1, 1, 1, 1, 1, 1, 2, 5),
                    byrow = T, nrow = 2)
design_weight <- design_ex[2,]/n
res_design_app <- tibble(x = design_app[1,], n_point = NA, weight = design_app[2,])
res_design_ex <- tibble(x = design_ex[1,], n_point = design_ex[2,], weight = design_weight)

calc_loss(res_design_app, 3) # 5.275251
calc_loss(res_design_ex, 3)  # 5.274603
# loss_appr <- 1.3188
# loss_ex <- 1.3187

# N 51 n18 ----------------------------------------------------------------
N <- 51; n <- 18; # seet 117
design_app <- matrix(c(-1, -0.44, 0.44, 1,
                       0.25, 0.25, 0.25, 0.25),
                     byrow=T,nrow=2)

design_ex <- matrix(c(-1, -0.4478, -0.4475, -0.4474, -0.4464, 0.4468, 0.4471, 0.4474,0.4475, 1,
                      5, 1, 2, 1, 1, 1, 1, 1, 1, 4),
                    byrow = T, nrow = 2)
design_weight <- design_ex[2,]/n
res_design_app <- tibble(x = design_app[1,], n_point = NA, weight = design_app[2,])
res_design_ex <- tibble(x = design_ex[1,], n_point = design_ex[2,], weight = design_weight)

calc_loss(res_design_app, 3) # 5.275251
calc_loss(res_design_ex, 3)  # 5.299447