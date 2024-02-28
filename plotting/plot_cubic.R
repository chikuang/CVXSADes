library(ggplot2)
library(dplyr)
library(tibble)
library(xtable)
library(latex2exp)
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 16),
             legend.text = element_text(size = 12),
             axis.text = element_text(size = 12),
             legend.title = element_blank(),
             title = element_text(size = 12))
loss <- matrix(c(1.3274, 1.3256,1.3224, 1.3188, 1.3188, 1.3187, 1.3187, 1.3187, 1.3187,
                   37.946,  38.915, 37.9, 37.589, 37.525, 37.525, 37.52, 37.52, 37.52,
                   5, 11, 21, 51, 101, 201, 501, 1001,2001), byrow = T, nrow = 3) %>% t()
colnames(loss) <- c("loss_D", "loss_A", "N")
loss_D_true <- 1.3187
loss_A_true <- 37.52

tbl_loss <- as_tibble(loss) %>% mutate(loss_D = exp(loss_D), 
                                       eff_D = exp(loss_D_true)/ loss_D,
                                       eff_A = loss_A_true/ loss_A) %>% 
  dplyr::select(N, loss_D, eff_D, loss_A, eff_A)

tbl_loss %>% 
  ggplot(aes(x = N, y = eff_D)) + 
  geom_line() +
  geom_point(size = 3) + 
  labs(title = TeX(" $Eff_D(\\xi_N^*)$  v.s. N"),
       x = "N", y = "D-efficiency") +
  ylim(0.95, 1)

tbl_loss %>% 
  ggplot(aes(x = N, y = eff_A)) + 
  geom_line() +
  geom_point(size = 3) + 
  labs(title = TeX(" $Eff_A(\\xi_N^*)$  v.s. N"),
       x = "N", y = "A-efficiency") +
  ylim(0.95, 1)
# print(xtable(tbl_loss_D, type = "latex"), file = "poly3_graph.tex")
xtable(tbl_loss, type = "latex",digits=c(0,0,4,4,4,4))


# Exact design ------------------------------------------------------------
loss <- matrix(c(1.3391, 1.3187, 1.3209, 1.3187,1.3195,
                 39.112, 38.193, 37.685, 37.693, 37.579,
                 10, 20, 30, 40,50),
               byrow = T, nrow = 3) %>% t()
colnames(loss) <- c("loss_D", "loss_A", "n")
loss_D_true <- 1.3187
loss_A_true <- 37.52
tbl_loss <- as_tibble(loss) %>% mutate(loss_D = exp(loss_D), 
                                       eff_D = exp(loss_D_true)/ loss_D,
                                       eff_A = loss_A_true/ loss_A) %>% 
  dplyr::select(n, loss_D, eff_D, loss_A, eff_A)

tbl_loss %>% 
  ggplot(aes(x = n, y = eff_D)) + 
  geom_line() +
  geom_point(size = 3) + 
  labs(title = TeX(" $Eff_D(\\xi_n^*)$  v.s. n"),
       x = "n", y = "D-efficiency") +
  ylim(0.95, 1)

tbl_loss %>% 
  ggplot(aes(x = n, y = eff_A)) + 
  geom_line() +
  geom_point(size = 3) + 
  labs(title = TeX(" $Eff_A(\\xi_n^*)$  v.s. n"),
       x = "n", y = "A-efficiency") +
  ylim(0.95, 1)
xtable(tbl_loss, type = "latex",digits=c(0,0,4,4,4,4)) 