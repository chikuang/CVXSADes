dat <- read.csv("./paper/plots/app3_c-opt_n12.csv",
                header = FALSE)
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

my_val <- unlist(dat)[-1] # the first one is the OAD
tibble(val = my_val) %>% rowid_to_column("iter") %>% 
  ggplot(aes(x = iter, y = val) ) + ylim(0, 0.25) + 
  scale_y_continuous(trans='log10') +
  geom_hline(yintercept = my_val[1], linetype="dashed", 
             color = "red", size=0.5) + 
  geom_line(col = "blue") + 
  labs(x = "iteration", 
       y = expression(paste("Loss function  ", Phi(xi[n]^(t)), "   (in log10 scale)"))) 
  

                     