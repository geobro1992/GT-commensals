##NMDS commensals#####

library(vegan)
library(tidyverse)
library(ggplot2)
library(lme4)

sites <- read.csv('sites_match.csv',
                  stringsAsFactors=FALSE)


commensals <- read.csv('Abundance_match.csv',
                       stringsAsFactors=FALSE,
                       na.strings=c('OUT','NA','UNK','',' ')) %>%
  tbl_df() %>%
  rename(Site.Name=Survey.Area) %>%
  left_join(., sites %>% 
              dplyr::select(Site.Name, Habitat.Type), by='Site.Name') %>%
  mutate_at(vars(Gopher.Frog:Turkey), 
            funs(ifelse(is.na(.), 0, .)))


com.dat = commensals %>%
  group_by(Site.Name, Burrow.ID, Habitat.Type) %>%
  summarise_at(vars(Gopher.Frog:Turkey),
               funs(sum(., na.rm=TRUE))) %>%
  ungroup()

ids = com.dat[,1:3]
counts = as.matrix(com.dat[,4:79])
counts[which(counts > 1)] = 1
y = rowSums(counts)

sum.dat = data.frame(ids, y)


group_ordered <- with(sum.dat,                       # Order boxes by median
                      reorder(Site.Name,
                              y,
                              median))

sum.dat$Site.Name = factor(sum.dat$Site.Name, levels = levels(group_ordered))


plot.dat = aggregate(y ~ Site.Name + Habitat.Type, data = sum.dat, function(x) quantile(x, probs = c(0.05, 0.5, 0.95)))
plot.dat = cbind(plot.dat[,1:2], plot.dat$y)
colnames(plot.dat) = c("Site", "Habitat", "lwr", "median", "upr")

g1 = ggplot(plot.dat, aes(x = Site, y = median)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr, color = Habitat)) +
  scale_color_manual(values=c('forestgreen','dark grey')) +
  ylab("Species Richness") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  theme_Publication()


ggsave("commensal_site_summary.png", plot = g1, width = 10, height = 5, units = "in", dpi = 600)



############
# NMDS
###########

commensals <- read.csv('Abundance_match.csv',
                       stringsAsFactors=FALSE,
                       na.strings=c('OUT','NA','UNK','',' ')) %>%
  tbl_df() %>%
  rename(Site.Name=Survey.Area) %>%
  left_join(., sites %>% 
              dplyr::select(Site.Name, Habitat.Type), by='Site.Name') %>%
  dplyr::select(-unk.warbler, -unk.frog, -unk.snake, -unk.lizard, -unk.mammal, 
                -unk.rodent, -unk.sparrow, -Canid, -unknown.bird, -unk.animal,
                -Burrowing.Owl) %>%
  mutate_at(vars(Gopher.Frog:Turkey), 
            funs(ifelse(is.na(.), 0, .)))



burrow.sizes <- commensals %>%
  mutate(width.numeric=as.numeric(str_extract(Width, "[0-9]+")),
         width.chr=str_extract(substr(Width, 1, 100), "[aA-zZ]+"),
         burrow.class1=ifelse(width.chr %in% c('L','adult','lg','large'), 'adult',
                              ifelse(width.chr %in% c('sm'), 'juvenile', 
                                     ifelse(width.chr %in% c('med','M','medium'), 
                                            'subadult', NA))),
         burrow.class2=ifelse(width.numeric<130, 'juvenile',
                              ifelse(width.numeric>=230, 'adult',
                                     ifelse(width.numeric>=130 &
                                              width.numeric<230,
                                            'subadult', NA))),
         burrow.class=as.factor(coalesce(burrow.class1, burrow.class2))) %>%
  dplyr::select(Habitat.Type, Site.Name, Season, Round, burrow.class) %>%
  group_by(Habitat.Type, Site.Name, Season, Round, burrow.class) %>%
  summarise(counts=n()) %>%
  ungroup() %>%
  na.omit() %>%
  spread(burrow.class, counts) %>%
  mutate(juvenile=ifelse(is.na(juvenile), 0, juvenile),
         subadult=ifelse(is.na(subadult), 0, subadult),
         total.burrows=adult + juvenile + subadult,
         prop.juvenile=juvenile/total.burrows)


site.spp <- commensals %>%
  group_by(Habitat.Type, Site.Name, Season, Round) %>%
  summarise_at(vars(Gopher.Frog:Turkey), funs(sum)) %>%
  ungroup() %>%
  left_join(., burrow.sizes,
            by=c('Habitat.Type','Site.Name','Season','Round'))

site.mds <- metaMDS(site.spp %>% 
                      dplyr::select(Gopher.Frog:Turkey),
                    distance='jaccard', binary=FALSE,
                    trymax=500)

treat=site.spp$Habitat.Type

png(filename = "GT_animal_NMDS.png", width = 7, height = 7, units = "in", res = 600)

par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
x <- c(-2, -1, 1, 2)
plot(x, c(-10, -10, -10, -10), type = "p", ylab = "", xlab = " ", cex = 1.5, 
     ylim = c(-1.5, 1.5), xlim = c(-1.5, 1.5), lwd = 2, pch = 5, axes = F, main = " ")
axis(1, at = c(-1, -0.5, 0, 0.5, 1.5))
mtext("NMDS1", side = 1, line = 3, cex = 1.5, font = 2)
axis(2, pos = 1.2)
par(las = 0)
mtext("NMDS2", side = 2, line = 2, cex = 1.5, font = 2)
ordiplot(site.mds,type="n", xlim = c(-1, 1), ylim = c(-1, 1.4))
ordiellipse(site.mds,groups=treat,draw="polygon",col=c('forestgreen','dark grey'),label=F,conf=0.95)
orditorp(site.mds,display="species",labels = FALSE, col="black",pch = 16, cex = 1)

#legend
points(0.15, 1.3, pch = 16, lwd = 2, cex = 2, col = "forestgreen")
text(0.3, 1.3, "Forested Site", cex = 1.5, font = 1, adj = 0)
points(0.15, 1.1, pch = 19, lwd = 2, cex = 2, col = "dark grey")
text(0.3, 1.1, "Test Range", cex = 1.5, font = 1, adj = 0)

dev.off()


################
# commensals GLM
################

sites <- read.csv('sites_match.csv',
                  stringsAsFactors=FALSE)


commensals <- read.csv('Abundance_match.csv',
                       stringsAsFactors=FALSE,
                       na.strings=c('OUT','NA','UNK','',' ')) %>%
  tbl_df() %>%
  rename(Site.Name=Survey.Area) %>%
  left_join(., sites %>% 
              dplyr::select(Site.Name, Habitat.Type), by='Site.Name') %>%
  mutate_at(vars(Gopher.Frog:Turkey), 
            funs(ifelse(is.na(.), 0, .)))

burrow.sizes <- commensals %>%
  mutate(width.numeric=as.numeric(str_extract(Width, "[0-9]+")),
         width.chr=str_extract(substr(Width, 1, 100), "[aA-zZ]+"),
         burrow.class1=ifelse(width.chr %in% c('L','adult','lg','large'), 'adult',
                              ifelse(width.chr %in% c('sm'), 'juvenile', 
                                     ifelse(width.chr %in% c('med','M','medium'), 
                                            'subadult', NA))),
         burrow.class2=ifelse(width.numeric<130, 'juvenile',
                              ifelse(width.numeric>=230, 'adult',
                                     ifelse(width.numeric>=130 &
                                              width.numeric<230,
                                            'subadult', NA))),
         burrow.class=as.factor(coalesce(burrow.class1, burrow.class2))) %>%
  dplyr::select(Habitat.Type, Site.Name, Season, Round, burrow.class) %>%
  group_by(Habitat.Type, Site.Name, Season, Round, burrow.class) %>%
  summarise(counts=n()) %>%
  ungroup() %>%
  na.omit() %>%
  spread(burrow.class, counts) %>%
  mutate(juvenile=ifelse(is.na(juvenile), 0, juvenile),
         subadult=ifelse(is.na(subadult), 0, subadult),
         total.burrows=adult + juvenile + subadult,
         prop.juvenile=juvenile/total.burrows)

com.dat <- commensals %>%
  group_by(Habitat.Type, Site.Name, Season, Round) %>%
  summarise_at(vars(Gopher.Frog:Turkey), funs(sum)) %>%
  ungroup() %>%
  left_join(., burrow.sizes,
            by=c('Habitat.Type','Site.Name','Season','Round'))


com.dat$alpha <- specnumber(com.dat[,5:80])
com.dat$shannon = diversity(com.dat[,5:80], index="shannon")

##################
# NOTE: RE gives singular fit (estimate = 0), so model result is the same when RE removed
m1 <- glmer(alpha ~ Season + Habitat.Type + total.burrows + prop.juvenile + (1|Site.Name), data = com.dat, family = "poisson")
summary(m1)
m1 <- glm(alpha ~ Season + Habitat.Type + total.burrows + prop.juvenile, data = com.dat, family = "poisson")
summary(m1)
library(report)
report(m1)

library(lmerTest)
m1 <- lmer(shannon ~ Season + Habitat.Type + total.burrows + prop.juvenile + (1|Site.Name), data = com.dat)
summary(m1)
report(m1)

library(nlme)
m1 <- lme(shannon~Season + Habitat.Type + total.burrows + prop.juvenile,random=~1|Site.Name,data=com.dat)
anova(m1)
