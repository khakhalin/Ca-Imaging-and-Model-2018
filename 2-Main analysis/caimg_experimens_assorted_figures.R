# ========================
# Reads data from a summary table (file names is hard-coded below), which in turn is exported
# from my summary Excel file. Then we can build different plots as we wish.
# ========================

require(tidyr)
require(dplyr)
require(ggplot2)

rm(list = ls())  # Clear workspace

# --- Network-level (brain-level) plots
d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/caimg_experiment_summary.csv",
                sep=',',header=T)
names(d)

# Encoding, between stages:
ggplot(d,aes(factor(Stage),encoding)) + theme_classic() + geom_point(shape=1)


# --- All-cell plots
d <- read.table("7_Ca imaging/git - CaImaging Paper/2-Main analysis/all_cells_latencies.csv",
                sep=',',header=T)
names(d)

ds = d %>% group_by(ibrain) %>% summarize(minlat = quantile(lat,0.1))
ggplot(d,aes(dist,lat)) + theme_classic() + geom_point(alpha=0.2)
# Now subtract quantiles from each lat, and add average (~250) back
# Then plut is as hex plot with alpha = density
