require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

d <- read.table("C:/Users/Arseny/Documents/3_Modeling/network_measures_for_r_180908.csv",sep=",",header=T)
names(d)

d2 <- gather(d,var,val,-Type,-Stage,-Name)
names(d2)

ggplot(subset(d2,var %in% c("Efficiency","Modul","Flow","cluster")),aes(Type,val,group=Name)) + theme_bw() +
  geom_line(color="lightblue") + geom_point(aes(color=Type)) + 
  facet_grid(var~Stage,scales="free_y")

# Assortativities look weird and a bit suspicious. Weird because after thorough rewiring they always
# converge to one value. Suspicious, because this value is not always 0 (IO and OI seem to converge to 0, 
# but II and OO converge to some non-zero negative value). I suspect that there's some bug-like 
# incompartibility between the weighted generalization of assortativity they use, and the way I randomize the 
# graphs. And as assortativity analysis seems to be a stretch and a dud overall, I am inclined to drop 
# it from the paper.
ggplot(subset(d2,var %in% c("II_assrt","OO_assrt","IO_assrt","OI_assrt")),aes(Type,val,group=Name)) + theme_bw() +
  geom_line(color="lightblue") + geom_point(aes(color=Type)) + 
  facet_grid(var~Stage,scales="free_y")
