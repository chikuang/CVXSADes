dat <- read.csv("./../app2-group/app3_c-opt_n12.csv")
library(readr)
library(ggplot2)
library(dplyr)
library(tibble)


theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5, size = 14),
             legend.text = element_text(size = 12),
             axis.text = element_text(size = 14),
             legend.title = element_blank(),
             title = element_text(size = 12))
tibble(val = unlist(dat)) %>% rowid_to_column("iter") %>% 
  ggplot(aes(x = iter, y = val) ) + ylim(0, 0.25) + 
  geom_line( col = "blue") + labs(x = "iteration", y = expression(paste("Loss function  ", Phi(xi[n]^(t)))))

                     