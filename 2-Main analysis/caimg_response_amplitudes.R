#
# Builds pictures for Fig 2, those about response amplitudes

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())      # Clear workspace

d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/avamps_allcells_allbrains.csv",
      sep=",",header=T)
names(d)
d$cell = 1:nrow(d)   # add cell ids
d2 = gather(d,type,amp,-ibrain,-stage,-cell)
head(d2)

ds = d2 %>% group_by(stage,ibrain,type) %>% 
  summarize(m = mean(amp), s = sd(amp))
ds$type = factor(ds$type,levels = c('f','s','c'))

ggplot(data=ds,aes(type,m*1000,group=ibrain,color=type)) + theme_classic() +
  geom_line(color='gray') + geom_point(shape=1) +
  facet_grid(.~stage)
