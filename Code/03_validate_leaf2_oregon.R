# Code files for EelgrassDiseaseTemperature manuscript
# 03_validate_leaf2_oregon

# Last updated 2022-03-30 by Lillian Aoki

# This script validates using leaf 2 for the Oregon measurements, instead of leaf 3 as in the other regions
# Outputs are Fig S8 in the manuscript
library(tidyverse)

or <- read.csv("Data/OR_compare_leaf2_leaf3.csv")
or$Estuary <- as.factor(or$Estuary)
or$Site <- as.factor(or$Site)
leaf <- tibble("Site"=or$Site,"Depth"=or$Depth,"Shoot"=or$`Shoot`,"Leaf"=or$`Leaf`,WD=or$`Wasting.Disease.`)
leaf <- na.omit(leaf)
leaf$WD <- gsub("1","0",leaf$WD)
leaf$WD <- gsub("2","1",leaf$WD)
leaf$Meadow[leaf$Site=="IF"] <- "OR_D"
leaf$Meadow[leaf$Site=="SBS"] <- "OR_C"
leaf$Meadow[leaf$Site=="SBN"] <- "OR_E"

leaf$Leaf <- as.factor(leaf$Leaf)
leaf_summ <- leaf %>%
  group_by(Site,Meadow,Depth,Leaf)%>%
  summarise(n=length(WD),diseased=length(WD[WD==1]),healthy=length(WD[WD==0]),
            prevalence=diseased/n)
leaf_summ2 <- leaf %>%
  group_by(Site,Meadow,Leaf)%>%
  summarise(n=length(WD),diseased=length(WD[WD==1]),healthy=length(WD[WD==0]),
            prevalence=diseased/n)
print(leaf_summ2)

ggplot(leaf,aes(x=Leaf,fill=WD))+geom_bar(position = position_dodge(preserve = "single"))+
  facet_wrap(~Meadow)+
  scale_y_continuous(expand = c(0,0),limits=c(0,62))+
  theme_bw(base_size = 11)+
  scale_fill_manual(values=c("darkgreen","grey50"),labels=c("Healthy","Diseased"))+
  xlab("Leaf rank")+
  ylab("Number of leaves")+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))
ggplot(leaf[leaf$Leaf==2|leaf$Leaf==3,],aes(x=Leaf,fill=WD))+geom_bar(position = position_dodge(preserve = "single"))+
  facet_wrap(~Meadow)+
  scale_y_continuous(expand = c(0,0),limits=c(0,62))+
  theme_bw(base_size = 11)+
  scale_fill_manual(values=c("darkgreen","grey50"),labels=c("Healthy","Diseased"))+
  xlab("Leaf rank")+
  ylab("Count of plants")+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white"))
ggsave(filename = "Figures/FigS1_leaf2_oregon.jpg", width = 4.75, height = 3.4)
# create high resolution version for submission (not uploaded)
ggsave(filename = "Figures/HighRes/FigS1_leaf2_oregon.tiff", width = 4.75, height = 3.4)

# check with logistic regression model
leaf23 <- subset(leaf,Leaf=="2" | Leaf=="3")
leaf23$WD <- as.integer(leaf23$WD)
m1 <- glm(WD~Leaf,data=leaf23, family=binomial)
summary(m1)
m2 <- glm(WD~Leaf*Site,data=leaf23, family=binomial)
summary(m2)
drop1(m2)
# no significant effect of Leaf 3