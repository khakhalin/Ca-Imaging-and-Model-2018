# ========================
# Script to analyze model outputs.
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

# Absolute address, replace as necessary
localFolder = 'C:/Users/Arseny/Documents/7_Ca imaging/git - CaImaging Paper/3-Model/'

### Combine several model outputs into one large dataframe
# 1st one creates the dataframe:
d <- read.table(paste(localFolder,"modelAnalysis181020 1 slide looming.csv",sep=''),sep=",",header=T)
d$type = 'Base'
# All others concatenate to it:
t <- read.table(paste(localFolder,"modelAnalysis181020 2 slide vis.csv",sep=''),sep=",",header=T)
t$type = 'Visual'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis181020 3 slide rand.csv",sep=''),sep=",",header=T)
t$type = 'Random'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis181023 4 slide looming nointrinsic.csv",sep=''),sep=",",header=T)
t$type = 'no Intrinsic'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis181023 5 slide looming Hebb.csv",sep=''),sep=",",header=T)
t$type = 'no STDP'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis181023 6 decay looming.csv",sep=''),sep=",",header=T)
t$type = 'no Competition'
d <- rbind(d,t)

names(d)
summary(d)

### Remove columns that were originally calculated, but that we don't like anymore:
d <- d %>% select (-c(nPCAto80,competition,sel90perc,sel90perc_SC,rCluSpk,rSelRnt,rSelGth,
                      shESelGrow,selEGrowth,nRichTo80F,clustPrefPval,revFlow,cycl,gammaIn,
                      rDirWei,rSelRnt,clusterCompactness)) 

dg <- gather(d,var,value,-file,-type,-stage,-rewire)
dg$var <- factor(dg$var)
levels(dg$var) # List of levels (different measures we have)

# Let's try to put measures in a somewhat meaningful order
levelSequence <- c(
  "file","type","stage",
  "fullBrainSel","meanSel","shareSelCells","sel90m50","bestPredict",
  "fullBrainSel_SC","meanSel_SC","shareSelCells_SC","sel90m50_SC",
  "rSelfcSelfs","rSelfcSelsc",
  "rPosSel","mDistWei",
  "rSelClu","rSelNet","rSelIns","selAssort","rSelSpk",
  "gammaIn","gammaOu","deg0","deg12","deg5p","recip",
  "nRichTo80C",
  "nClusters","clusterPreference",
  "eff","modul","clust","flow")
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

# Plot Averages only, many TYPES in each plot
ggplot(dgs,aes(stage,m,color=type)) + 
  geom_point(size=1,shape=1) + 
  geom_line(aes(group=type)) +
  facet_wrap(~var,scales = "free_y") +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

dgs %>% filter(stage==5,var=="recip")

# ----------- end of meaningful analyses -------------

# Everything below are old attempts when I tried to compare "real" values to values on a rewired
# graph. It's really hard to visualize, and even harder to interpret, so these attempts are
# abandoned, at least for now.


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
