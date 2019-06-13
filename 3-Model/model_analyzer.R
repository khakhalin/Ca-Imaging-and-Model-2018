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
d <- read.table(paste(localFolder,"modelAnalysis 1 190322 slide looming.csv",sep=''),sep=",",header=T)
d$type = 'Base'
# All others concatenate to it:
t <- read.table(paste(localFolder,"modelAnalysis 2 190322 slide vis.csv",sep=''),sep=",",header=T)
t$type = 'Visual'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis 3 190322 slide rand.csv",sep=''),sep=",",header=T)
t$type = 'Random'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis 4 190322 slide looming nointrinsic.csv",sep=''),sep=",",header=T)
t$type = 'no Intrinsic'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis 5 190611 slide looming Hebb.csv",sep=''),sep=",",header=T)
t$type = 'no STDP'
d <- rbind(d,t)
t <- read.table(paste(localFolder,"modelAnalysis 6 190322 decay looming.csv",sep=''),sep=",",header=T)
t$type = 'no Competition'
d <- rbind(d,t)

names(d)
summary(d)

nExp <- nrow(d)/5
d$exp <- rep(seq(1,nExp),each=5) # Label experiments

# Remove columns that are still calculated by the scripts, but that we don't like anymore:
d <- d %>% dplyr::select (-c(nPCAto80,competition,sel90perc,sel90perc_SC,rCluSpk,rSelRnt,rSelGth,
                      shESelGrow,selEGrowth,nRichTo80F,clustPrefPval,revFlow,cycl,gammaIn,
                      rDirWei,rSelfcSelfs,rSelfcSelsc,deg0,deg12,deg5p,
                      clusterPreference,clusterCompactness)) 

# -- Correct weird values
# For some reason for "no competition" the starting point (the very first one) is weird.
# something is wrong with renormalization of synapses for this very first snapshot.
# So tucking these values in.
d$clust <- pmin(d$clust,0.02) 
ggplot(d,aes(stage,clust,color=type,group=file)) + theme_minimal() + geom_line()

d$eff[d$stage==1] <- pmin(d$eff[d$stage==1],0.04)
ggplot(d,aes(stage,eff,color=type,group=file)) + theme_minimal() + geom_line()

d$gammaOu <- abs(d$gammaOu) # from negative to positive values


# ---------- Gather, to prepare for summary statistics bellow
dg <- gather(d,var,value,-file,-type,-stage,-rewire,-exp)
dg$var <- factor(dg$var)
levels(dg$var) # List of levels (different measures we have)

# Let's try to put measures in a somewhat meaningful order
levelSequence <- c(
  "file","type","stage",
  "fullBrainSel","meanSel","shareSelCells","sel90m50","bestPredict",
  "fullBrainSel_SC","meanSel_SC","shareSelCells_SC","sel90m50_SC",
  "rSelfcSelfs","rSelfcSelsc",
  "rPosSel","mDistWei",
  "rSelSpk","rSelClu","rSelNet","rSelIns","selAssort",
  "nRichTo80C","nClusters","gammaOu","recip",
  "synfire","synHelp",
  "clusterPreference","clusterCompactness",
  "eff","modul","clust","flow")
existingLevels <- levels(dg$var)
dg$var <- factor(dg$var,levels=intersect(levelSequence,existingLevels))
levels(dg$var)

# Summary data frame
dgs = dg %>% group_by(type,stage,var,rewire) %>% summarize(
  m = mean(value),
  n = n(),
  s = sd(value),
  ci = -s/sqrt(n)*qt(0.025,df=n-1)) # Averages and cis
head(dgs)

# Plot Averages only, several MODEL TYPES in each plot <---- Main plot
ggplot(dgs,aes(stage,m,color=type,linetype=type,size=type)) + 
  geom_point(size=1,shape=1) + 
  geom_line(aes(group=type)) +
  facet_wrap(~var,scales = "free_y", ncol=5) +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_linetype_manual(values=c("solid","solid","solid","solid","dashed","dashed"))+
  scale_color_manual(values=c('black','red','palegreen3','dodgerblue','violetred','gray'))+
  scale_size_manual(values=c(1, 0.5, 0.5, 0.5, 0.5, 0.5))+
  NULL

# Values by the end of simulation
names(dgs)
levels(dgs$var)
dgs %>% dplyr::filter(stage==5,var=="meanSel")
dgs %>% dplyr::filter(stage==5,var=="shareSelCells")
dgs %>% dplyr::filter(stage==5,var=="meanSel_SC")
dgs %>% dplyr::filter(stage==5,var=="shareSelCells_SC")
dgs %>% dplyr::filter(stage==5,var=="rSelfcSelsc")

# Same, but only global clustering
ggplot(subset(dgs,var=="clust"),aes(stage,m,color=type)) + 
  geom_point(size=1,shape=1) + 
  geom_line(aes(group=type)) +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# Same, but only global efficiency
ggplot(subset(dgs,var=="eff"),aes(stage,m,color=type)) + 
  geom_point(size=1,shape=1) + 
  geom_line(aes(group=type)) +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL


# --- Synfire analysis
# All curves for one type of experiments, and one var
ggplot(subset(dg,type=="Base" & var %in% c("synfire")),aes(stage,value,group=exp)) + 
  geom_point(size=1,shape=1) + 
  geom_line() +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# Across types
ggplot(subset(dgs,var %in% c("synfire","synHelp")),aes(stage,m,color=type)) + 
  geom_point(size=1,shape=1) + 
  geom_line(aes(group=type)) +
  theme_bw() +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~var,scales = "free_y") +
  NULL

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
