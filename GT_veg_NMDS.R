##GT NMDS veg#####

library(vegan)
library(tidyverse)

veg<-read.csv('herbaceousFull.csv',header = T,row.names = 1)

m<-metaMDS(veg,k=4,distance='bray',trymax=100)
stressplot(m)
treat=c(rep("Forested",7),rep("Test Range",10))

png(filename = "GT_veg_NMDS.png", width = 7, height = 7, units = "in", res = 600)

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
x <- c(1, 2, 3, 4)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = "", xlab = " ", cex = 1.5, 
     ylim = c(0.3, 0.6), xlim = c(1, 4), lwd = 2, pch = 5, axes = F, main = " ")
axis(1, at = c(1.5, 2.5, 3.5), labels = c("HF", "LF", "VLF"))
mtext("NMDS1", side = 1, line = 3, cex = 1.5, font = 2)
axis(2, pos = 1.2)
par(las = 0)
mtext("NMDS2", side = 2, line = 2, cex = 1.5, font = 2)
ordiplot(m,type="n")
ordiellipse(m,groups=treat,draw="polygon",col=c('forestgreen','dark grey'),label=F,conf=0.95)
orditorp(m,display="species",labels = FALSE, col="black",pch = 16, cex = 1)

#legend
points(0.65, 1.3, pch = 16, lwd = 2, cex = 2, col = "forestgreen")
text(0.8, 1.3, "Forested Site", cex = 1.5, font = 1, adj = 0)
points(0.65, 1.1, pch = 19, lwd = 2, cex = 2, col = "dark grey")
text(0.8, 1.1, "Test Range", cex = 1.5, font = 1, adj = 0)

dev.off()
