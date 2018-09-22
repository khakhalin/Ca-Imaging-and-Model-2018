# Script to analyze model outputs

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

d <- read.table("C:/Users/Arseny/Documents/3_Modeling/modelAnalysis180918 slide looming.csv",sep=",",header=T)

names(d)

d <- d %>% select (-c(asII,asIO,asOI,asOO,rDistWei,nPCAto80)) # Drop assortativities as we don't like them anymore
      # And some other useless measures as welll

dg <- gather(d,var,value,-file,-type,-competition,-stage)
# Put measures in a somewhat meaningful order:
dg$var <- factor(dg$var,levels=c(
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
        "eff","modul","clust","flow","revFlow","cycl"))
head(dg)

dgs = dg %>% group_by(type,stage,var) %>% summarize(
  m = mean(value),
  n = n(),
  s = sd(value),
  ci = -s/sqrt(n)*qt(0.025,df=n-1)) # Averages and cis
head(dgs)

# Averages
ggplot(dgs,aes(stage,m,color=type)) + theme_bw() +
  geom_point() + geom_line(aes(group=type)) +
  facet_wrap(~var,scales = "free_y") +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# All values and averages
ggplot(data=dg,aes(stage,value,color=type)) + theme_bw() +
  geom_point(alpha=0.5) + geom_line(aes(group=file),alpha=0.2) +
  facet_wrap(~var,scales = "free_y") +
  geom_point(data=dgs,aes(stage,m),color="black") +
  geom_line(data=dgs,aes(stage,m,group=type),color="black") +
  theme(axis.text.y=element_text(size=6),
        strip.background=element_rect(linetype='blank',fill='white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  NULL

# Zoom in on cyclicity in particular
ggplot(data=d,aes(stage,cycl,color=type)) + theme_bw() +
  geom_point() + geom_line(aes(group=file)) +
  scale_y_log10() +
  NULL

names(dgs)
dgs2 <- subset(dgs,stage==5)
data.frame(sprintf('%s - %4.2f pm %4.2f',dgs2$var,dgs2$m,dgs2$s)) # means and sd by the end of training
