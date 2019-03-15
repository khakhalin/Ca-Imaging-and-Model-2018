# ========================
# Reads network measures (file names is hard-coded below), for every brain, both actual,
# and calculated on rewired (scrambled) graphs. Plots them.
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

# --------------------------------------------------------------------------
# -------------- Network measures (with and without rewiring) --------------
# --------------------------------------------------------------------------
# For other experiment-level summaries that do not need rewiring, see below.

d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/caimg_network_measures.csv",sep=",",header=T)
names(d)

d2 <- gather(d,var,val,-type,-stage,-name)
names(d2)

ggplot(subset(d2,var %in% c("eff","modul","flow","clust")),aes(type,val,group=name)) + 
theme_bw() + 
 theme(panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
       axis.line = element_line()) + 
  geom_line(color="gray90") + 
  geom_point(aes(color=stage),shape=1) + 
  facet_grid(var~stage,scales="free_y") +
  stat_summary(aes(group=type),fun.y="mean",geom="point",shape=0,color="black") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  NULL

# Assortativities look weird and a bit suspicious. Weird because after thorough rewiring they always
# converge to one value. Suspicious, because this value is not always 0 (IO and OI seem to converge to 0, 
# but II and OO converge to some non-zero negative value). I suspect that there's some bug-like 
# incompartibility between the weighted generalization of assortativity they use, and the way I randomize the 
# graphs. And as assortativity analysis seems to be a stretch and a dud overall, I am inclined to drop 
# it from the paper.
ggplot(subset(d2,var %in% c("II_assrt","OO_assrt","IO_assrt","OI_assrt")),aes(Type,val,group=Name)) + theme_bw() +
  geom_line(color="lightblue") + geom_point(aes(color=Type)) + 
  facet_grid(var~Stage,scales="free_y")


# -----------------------------------------------------------
# -------------- Simple experiment statistics (different data file)

rm(list = ls())  # Clear workspace

dbasic <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/caimg_experiment_summary.csv",sep=",",header=T)
names(dbasic)

dbasic2 <- gather(dbasic,var,val,-stage,-name)
names(dbasic2)

# Degree distribution power:
ggplot(dbasic,aes(factor(stage),-gamma)) + 
  theme_classic() + 
  geom_point(aes(color=factor(stage)),shape=1,size=3) + 
  stat_summary(fun.y="mean",geom="point",shape=0,color="black") +
  NULL

# Assortativity of selectivity:
ggplot(data=dbasic,aes(factor(stage),selassort)) + 
  theme_classic() +
  geom_jitter(h=0,w=0.1,shape=1,size=2.5) + 
  NULL

# Several measurements related to clusters
ggplot(data=subset(dbasic2,var %in% c("nEns","maxModul","clustCompact")),
       aes(factor(stage),val)) + 
  theme_classic() +
  geom_jitter(h=0,w=0.1,shape=1,size=2.5) + 
  facet_wrap(~var,scales="free") +
  NULL
t.test(data=dbasic,nEns~stage) # 0.9
t.test(data=dbasic,maxModul~stage) # 0.03
t.test(data=dbasic,clustCompact~stage) # 0.3
