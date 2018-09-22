# Various analyses for the Ca imaging 2018 paper

# ------- Analysis of positions
# (Excel table "Selectivity analysis", tab "pos_full")

d = read.table("clipboard",header=T)

d$stage = factor(d$stage)
ds <- d %>% group_by(stage,name) %>% summarize(m=mean(sel))
summary(aov(data=d,sel~stage+name))   # On all points there is a difference
summary(aov(data=ds,m~stage))         # But no difference on summary points.

# Which one to believe? Let's visualize all points:
ggplot(d,aes(sel,name,color=stage)) + theme_bw() + geom_point() + 
  facet_grid(stage~.,scales='free_y') + geom_point(data=ds,aes(m,name),color='red')
# And summaries alone:
ggplot(ds,aes(factor(stage),m)) + theme_bw() + geom_point()
# See, everything hinges on these two really weird experiments. And one of them (with 
# mostly-negative selectivities) is particularly weird. So we cannot believe full data,
# and should stick with non-significant summary test.

# But what we can do is cancel effect of individual experiments for position analysis:
summary(aov(data=d,sel~stage+name+toOrigin*stage))
# Now visualize:
fit = lm(data=d,sel~name)
d$res = resid(fit)
ggplot(d,aes(toOrigin,res,color=stage)) + geom_point(alpha=0.3) + theme_bw() +
  geom_smooth(method=lm,se=F,color='black') + facet_wrap(~stage) +
  ylab('Adjusted selectivity')

# Lateral position effect is significant even after compensations:
summary(aov(data=d,sel~stage+name+toOrigin+x*stage))
fit = lm(data=d,sel~name+toOrigin)
d$res = resid(fit)
ggplot(d,aes(x,res,color=stage)) + geom_point(alpha=0.3) + theme_bw() +
  geom_smooth(method=lm,se=F,color='black') + facet_wrap(~stage) +
  ylab('Adjusted selectivity')
# Correlations for individual experiments
ds <- d %>% group_by(stage,name) %>% summarize(rho=cor(x,res),
                                               pval=cor.test(x,res)$p.value)
# How many were significant
ds %>% group_by(stage) %>% summarize(n=sum(pval<0.05))
# Unfortunately they were significant in both directions, so it makes no sense:
ggplot(ds,aes(rho,log(pval))) + geom_point() + theme_bw()
# Conclusion: probably an artifact, and I definitely don't believe in total result.
