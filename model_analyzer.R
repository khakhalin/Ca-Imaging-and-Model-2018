# ========================
# Script to analyze model outputs.
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

#d <- read.table("C:/Users/Arseny/Documents/7_Ca imaging/Model analysis/modelAnalysis180919 slide looming - with 50 rewirings for each.csv",sep=",",header=T)

### Combine several model outputs into one large dataframe
d <- read.table("C:/Users/Arseny/Documents/7_Ca imaging/Model analysis/modelAnalysis180925 1 slide looming updated measures.csv",sep=",",header=T)
t <- read.table("C:/Users/Arseny/Documents/7_Ca imaging/Model analysis/modelAnalysis180929 2 slide vis.csv",sep=",",header=T)
d <- rbind(d,t)
t <- read.table("C:/Users/Arseny/Documents/7_Ca imaging/Model analysis/modelAnalysis180929 3 slide shuffle.csv",sep=",",header=T)
d <- rbind(d,t)
t <- read.table("C:/Users/Arseny/Documents/7_Ca imaging/Model analysis/modelAnalysis180929 4 slide random.csv",sep=",",header=T)
d <- rbind(d,t)

names(d)
summary(d)

### Remove columns that we calculated, but that we don't like anymore:
# For old files: drop assortativities (useless), as well as some other low-impact measures
# d <- d %>% select (-c(asII,asIO,asOI,asOO,rDistWei,nPCAto80)) 
d <- d %>% select (-c(nPCAto80)) 

dg <- gather(d,var,value,-file,-type,-competition,-stage,-rewire)
dg$var <- factor(dg$var)
levels(dg$var) # List of levels (different measures we have)

# Let's try to put measures in a somewhat meaningful order
levelSequence <- c(
  "file","type","competition","stage",
  "fullBrainSel","meanSel","shareSelCells","sel90perc","sel90m50",
  "sumOTpredict","bestPredict",
  "fullBrainSel_SC","meanSel_SC","shareSelCells_SC","sel90perc_SC","sel90m50_SC",
  "rPosSel","rPosInf","rSelInf",
  "rDirWei","rDistWei","mDistWei",
  "rSelClu","rCluSpk","rSelNet","rSelRnt",
  "rSelGth","rSelIns","selAssort","shESelGrow","selEGrowth",
  "gammaIn","gammaOu","deg0","deg12","deg5p",
  "nRichTo80F","nRichTo80C","nPCAto80",
  "nClusters","clustVarExplained","clusterPreference","clusterCompactness","clustPrefPval",
  "rSelSpk","rSelfcSelfs","rSelfcSelsc",
  "eff","modul","clust","flow","revFlow","cycl")
existingLevels <- levels(dg$var)
dg$var <- factor(dg$var,levels=intersect(levelSequence,existingLevels))
head(dg)

# Summary data frame
dgs = dg %>% group_by(type,stage,var,rewire) %>% summarize(
  m = mean(value),
  n = n(),
  s = sd(value),
  ci = -s/sqrt(n)*qt(0.025,df=n-1)) # Averages and cis
head(dgs)

# Plot Averages only, many TYPES, but no rewire
ggplot(dgs,aes(stage,m,color=type)) + theme_bw() +
  geom_point() + geom_line(aes(group=type)) +
  facet_wrap(~var,scales = "free_y") +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# Plot Averages only, one type, but with REWIRE
ggplot(dgs,aes(stage,m,color=rewire)) + theme_bw() +
  geom_point() + geom_line(aes(group=rewire)) +
  facet_wrap(~var,scales = "free_y") +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# All values and averages for one TYPE, but with rewires. Take a while to render
ggplot(data=dg,aes(stage,value,color=rewire)) + theme_bw() +
  geom_point(alpha=0.5) + 
  geom_line(data=subset(dg,rewire=="original"),aes(group=file),alpha=0.3) +
  facet_wrap(~var,scales = "free_y") +
  geom_point(data=dgs,aes(stage,m),color="black") +
  geom_line(data=dgs,aes(stage,m,group=rewire),color="black") +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL
# It's impossible to connect rewired points; rewires are random, so we have to subset geom_line(),
# to only connect "proper" (original) points.

# Zoom in on cyclicity in particular
ggplot(data=d,aes(stage,cycl,color=rewire)) + theme_bw() +
  geom_point(alpha=0.5) + 
  geom_line(data=subset(d,rewire=="original"),aes(group=file),alpha=0.3) +
  scale_y_log10() +
  NULL

names(dgs)
dgs2 <- subset(dgs,stage==5)
data.frame(sprintf('%s - %4.2f pm %4.2f',dgs2$var,dgs2$m,dgs2$s)) # means and sd by the end of training
