# make a plot for the T-shirts/coffee cups


library(ggplot2)
library(ordbetareg)
library(tidyverse)
library(ggdist)
library(distributional)
library(ggfortify)

ggdistribution(dordbeta, seq(0, 1, 0.1), mu=0.5,phi=0.7,
               cutpoints=c(-2,2))


y <- rordbeta(10000,mu=0.5,phi=1,cutpoints=c(-2,2))


tibble(x=y) %>%
  ggplot(aes(x=x)) +
  #stat_slab(aes(x), density = "bounded", fill = "gray", color = "black") +
  geom_histogram(fill="white") +
  theme_ggdist() +
  labs(y="",x="") +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.text.x=element_text(face="bold",size=35,colour="white"),
        axis.line = element_line(colour="#050a30",linewidth=0),
        panel.border=element_rect(fill=NA,colour=NA),
        panel.background = element_rect(fill=NA,colour=NA),
        plot.background = element_rect(fill=NA,colour=NA))

# "#050a30"

  ggsave("tshirt_coffe.png",width=12,height=10,dpi=1000)
