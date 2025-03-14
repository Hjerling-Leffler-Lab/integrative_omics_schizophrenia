library(data.table)
library(ggplot2)
library(patchwork)

setwd("./9_sc_eQTL/output/fig4F/")

dat <- fread("df_3CT0.05_vs_GTEx.tsv.gz") 
# Note on the data:
# The eQTL pairs were mapped across the three cell classes. An eQTL significant in any cell class was kept in the dataset.

p.ei <- ggplot(dat, 
                aes(x=Excitatory, y=Inhibitory) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat$Inhibitory, dat$Excitatory),3),sep=""))
p.en <- ggplot(dat, 
                aes(x=Excitatory, y=NonNeuronal) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat$Excitatory, dat$NonNeuronal),3),sep=""))
p.in <- ggplot(dat, 
                aes(x=Inhibitory, y=NonNeuronal) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  xlim(axislow,axishigh) + ylim(axislow,axishigh) +
  labs(title=paste("corr=",round(cor(dat$Inhibitory, dat$NonNeuronal),3),sep=""))

p <- (p.ei + guides(fill="none")) | (p.en + guides(fill="none")) |(p.in + guides(fill="none"))
ggsave(plot=p,
       filename=paste(outdir,"fig4F_corr_btw_3CT.pdf",sep=""),
       width=9, height = 3)
