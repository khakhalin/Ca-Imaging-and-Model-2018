# ========================
# Builds pictures for Caimg paper Fig 2, about response selectivity
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())      # Clear workspace

d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/sel_allcells_allbrains stable.csv",
      sep=",",header=T)
names(d)
d$cell = 1:nrow(d)   # add cell ids
names(d)
d2 = gather(d,type,amp,-ibrain,-stage,-cell)
head(d2)

ds <- d2 %>% group_by(stage,ibrain,type) %>% 
  summarize(m = mean(amp), s = sd(amp))                # Average selectivity for each brain
ds$type = factor(ds$type,levels = c('fc','fs','sc'))   # Correct sequence, for the plot
ds <- subset(ds,type %in% c('fc','sc'))                # drop fs, we don't report it

# Average selectivity for each brain. Isn't used right now.
ggplot(data=ds,aes(type,m*1000,group=ibrain,color=type)) + theme_classic() +
  geom_line(color='gray') + geom_point(shape=1) +
  facet_grid(.~stage)
