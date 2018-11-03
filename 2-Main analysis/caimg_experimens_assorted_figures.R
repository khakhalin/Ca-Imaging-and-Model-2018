# ========================
# Reads data from a summary table (file names is hard-coded below), which in turn is exported
# from my summary Excel file. Then we can build different plots as we wish.
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/caimg_network_summary.txt",
                sep='\t',header=T)
names(d)

ggplot(d,aes(factor(Stage),encoding)) + theme_classic() + geom_point(shape=1)
