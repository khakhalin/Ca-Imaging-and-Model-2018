# ========================
# Caimg Centrality Analysis
# ========================

# Reads sel_centrality_allcells.csv, makes the plots

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/sel_centrality_allcells.csv",sep=",",header=T)
names(d)

# --- In-degrees ---
# All plots in facets
ggplot(data=d,aes(indegree,sel)) + 
  theme_bw() +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm",se=F,size=0.1) +
  facet_wrap(~nexp) +
  NULL
# Strong correlations: 1 2 5 7 16

# Selected plot to llustrate the correlation:
ggplot(data=subset(d,nexp==16),aes(indegree,sel)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL

# Summary plots
ds <- d %>% group_by(nexp) %>% summarise(cis = cor(sel,indegree))
ds$stage = ifelse(ds$nexp<15,46,49) # First 14 experiments are from younger tadpoles

ds

# THIS NEEDS TO BE FIXED
# Quick summary for now: values here do not match those in Excel (tab rSelNet)
# And those values from Excel are not saved in any csv for now.
# (They are at experiment level, so may be saved in network_measures file)
# But first we need to understand why the difference, by comparing Matlab code to this code.

t.test(data=ds,cis~stage) # No change in development, p=0.5
t.test(subset(ds,stage==46)$cis,mu=0)
t.test(subset(ds,stage==49)$cis,mu=0)

ggplot(data=ds,aes())
