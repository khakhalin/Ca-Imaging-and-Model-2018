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

# --- Construct summary dataset
ds <- d %>% group_by(nexp) %>% summarise(
  cis = cor(selFC,indegree),
  cks = cor(selFC,katz),
  css = cor(selFC,spiking),
  csss = cor(selFC,spiking,method="spearman")
  )
ds$stage = ifelse(ds$nexp<15,46,49) # First 14 experiments are from younger tadpoles

ds


# -------- In-degree --------
# All plots in facets
ggplot(data=d,aes(indegree,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm",se=F,size=0.1) +
  facet_wrap(~nexp) +
  NULL
# Strong correlations: 1 2 5 7 8 16 23

# Selected plot to llustrate the correlation:
ggplot(data=subset(d,nexp==16),aes(indegree,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL

# Al points ever observed
ggplot(data=d,aes(indegree,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL
cor.test(data=d,~selFC+indegree) # 0.3 - no effect for this extreme overpower haha

# Comparisons:
t.test(data=ds,cis~stage) # No change in development, p=0.5
t.test(subset(ds,stage==46)$cis,mu=0) # 0.09
t.test(subset(ds,stage==49)$cis,mu=0) # 0.14
t.test(ds$cis,mu=0) # 0.02

ggplot(data=ds,aes(factor(stage),cis)) + theme_bw() +
  geom_jitter(h=0,w=0.1,shape=1,size=2)



# ----------- Katz centrality -----------
# Selected plot to llustrate the correlation:
ggplot(data=subset(d,nexp==16),aes(katz,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL

# All points
ggplot(data=d,aes(katz,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL
cor.test(data=d,~selFC+katz) # 1e-6

# All points for a raster figure
ggplot(data=d,aes(katz,selFC)) + 
  theme_classic() +
  geom_point(alpha=0.1) +
  NULL

# Same, but with a fancy scatterplot (not clear enough for the paper though)
ggplot(data=d,aes(katz,selFC)) + 
  theme_bw() +
  stat_binhex(aes(alpha=..count..),fill="black") +
  geom_smooth(method="lm",se=F,size=0.1) +
  NULL

# Katz correlates with in-degree, but in a very weird way: two typical slopes
# Not related to age. Some sort of really unfortunate data variability, at topological level.
ggplot(data=d,aes(indegree,katz,color=factor(nexp<15))) +
  theme_bw() + geom_point(alpha=0.2) +
  NULL

# Comparisons:
t.test(data=ds,cks~stage) # No change in development, p=0.5
t.test(subset(ds,stage==46)$cks,mu=0) # 0.09
t.test(subset(ds,stage==49)$cks,mu=0) # 0.16
t.test(ds$cks,mu=0) # 0.03

ggplot(data=ds,aes(factor(stage),cks)) + 
  theme_bw() +
  geom_jitter(h=0,w=0.1,shape=1)


# ----- Spiking -----

# Total. Weird outlier for s46 (experiment #9)
ggplot(data=ds,aes(factor(stage),css)) + 
  theme_bw() +
  geom_jitter(h=0,w=0.1,shape=1)

# All plots in facets
ggplot(data=d,aes(spiking,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.1) +
  geom_smooth(method="lm",se=F,size=0.1) +
  facet_wrap(~nexp) +
  NULL
# Outliers in 2 experiments (one up, one down); maybe even 4

# Weird experiment 9 has a negative correlation because of one point only!
ggplot(data=subset(d,nexp==9),aes(spiking,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  geom_smooth(method="lm",se=F,size=0.1)

# We can remove this outlier:
ggplot(data=subset(d,nexp==9 & spiking<0.02),aes(spiking,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.8) +
  geom_smooth(method="lm",se=F,size=0.1)

# All points together:
ggplot(data=d,aes(spiking,selFC)) + 
  theme_bw() +
  geom_point(alpha=0.2) +
  geom_smooth(method="lm",se=F,size=0.1)
cor.test(data=d,~selFC+spiking) # p=0.006

t.test(data=ds,css~stage) # Almost a change in development, p=0.06
t.test(data=subset(ds,nexp!=9),css~stage) # With this one omitted, p=0.003
t.test(subset(ds,stage==46)$css,mu=0) # 0.001
t.test(subset(ds,stage==49)$css,mu=0) # 0.07
t.test(ds$css,mu=0) # 0.0003

# - Let's switch to Spearman correlations:
cor.test(data=d,~selFC+spiking,method="spearman") # p=4e-4
t.test(data=ds,csss~stage) # p=0.02
ds %>% group_by(stage) %>% summarize(m=mean(csss),s=sd(csss))
t.test(subset(ds,stage==46)$csss,mu=0) # 2e-5
t.test(subset(ds,stage==49)$csss,mu=0) # 0.04

ggplot(data=ds,aes(factor(stage),csss)) + 
  theme_classic() +
  geom_jitter(h=0,w=0.1,shape=1,size=2.1)



# -------- Local Clustering ---------

cor.test(data=d,~clust+selFC) # 0.6
nrow(d)
