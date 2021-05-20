# Code files for EelgrassDiseaesTemperature manuscript
# 01_eelgrass_survey_metrics

# Last updated 2021-05-20 by Lillian Aoki

# This script imports and visualizes data from the eelgrass wasting disease surveys in 2019
# Outputs include Fig 6 and Fig S5 in the manuscript

#### Load libraries ####
library(tidyverse)
library(plotrix)
library(cowplot)
region_order <- c("AK","BC","WA","OR","BB","SD")

# Read in data ####
## site coordinates
coords <- read.csv("Data/site_lat_long.csv")
# summarize by region, subtract 5 from mean latitude so that the region mean will plot to the right of all sites
coords_region <- coords %>%
  group_by(Region) %>%
  summarise(Lat=mean(Lat)-5,Long=mean(Long))
# in figures, the data are arranged by latitude to show the geographic gradient 

# survey metrics at the site level
all <- read.csv("Data/all_survey_metrics_site.csv")
all$Type <- "Site"
all <- left_join(all, coords, by=c("Region","SiteCode"))
# calculate canopy height as the combination of sheath length and longest blade length
all$CanopyHeightMean <- all$LongestBladeLengthMean+all$SheathLengthMean
all$CanopyHeightSe <- sqrt((all$LongestBladeLengthSe)^2+(all$SheathLengthSe)^2)

all_region <- all %>%
  group_by(Region, Type="Region")%>%
  summarise(DensityShootsSe=sd(DensityShootsMean)/sqrt(length(DensityShootsMean)),
            DensityShootsMean=mean(DensityShootsMean,na.rm=TRUE),
            EpiphytePerAreaSe=sd(EpiphytePerAreaMean)/sqrt(length(EpiphytePerAreaMean)),
            EpiphytePerAreaMean=mean(EpiphytePerAreaMean),
            LongestBladeLengthSe=sd(LongestBladeLengthMean)/sqrt(length(LongestBladeLengthMean)),
            LongestBladeLengthMean=mean(LongestBladeLengthMean),
            SheathLengthSe=sd(SheathLengthMean)/sqrt(length(SheathLengthMean)),
            SheathLengthMean=mean(SheathLengthMean),
            CanopyHeightSe=sd(CanopyHeightMean)/sqrt(length(CanopyHeightMean)),
            CanopyHeightMean=mean(CanopyHeightMean),
            PrevalenceSe=sd(PrevalenceMean)/sqrt(length(PrevalenceMean)),
            PrevalenceMean=mean(PrevalenceMean),
            SeveritySe=sd(SeverityMean)/sqrt(length(SeverityMean)),
            SeverityMean=mean(SeverityMean),
            BladeAreaSe=sd(BladeAreaMean)/sqrt(length(BladeAreaMean)),
            BladeAreaMean=mean(BladeAreaMean))

all_region <- left_join(all_region,coords_region,by=c("Region"))

all_region$SiteCode <- "Rg"
all_region$Region <- ordered(all_region$Region,levels=region_order)
all$Region <- ordered(all$Region,levels=region_order)
all <- dplyr::select(all,c("SiteCode",names(all_region)))

# combine site and region data for plotting
combo <- rbind(all,all_region)
combo$Region <- ordered(combo$Region, levels=region_order)
combo <- combo %>%
  ungroup%>%
  arrange(Region,-Lat)%>%
  mutate(order=row_number())

# Create plots ####

# shoot density plot
p_den <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),y=DensityShootsMean,ymin=DensityShootsMean-DensityShootsSe,
                                 ymax=DensityShootsMean+DensityShootsSe,color=Type,shape=Type))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","green4"))+
  scale_shape_manual(values=c(8,16))+
  scale_y_log10(limits=c(10,3000))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab(expression(atop("Shoot density",paste("(shoots m"^"-2"~")"))))+
  xlab("Site within Region")+
  theme_bw()+
  #theme(axis.text.x=element_text(angle=60,vjust = -0.05,hjust=0.25))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        legend.position = "none")
p_den

# epiphyte load plot
p_epi <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),y=EpiphytePerAreaMean,
                                     ymax=EpiphytePerAreaMean+EpiphytePerAreaSe,
                                     ymin=EpiphytePerAreaMean-EpiphytePerAreaSe,
                                     color=Type,shape=Type))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","green4"))+
  scale_shape_manual(values=c(8,16))+
  #scale_y_log10(limits=c(0.0000001,0.1))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab(expression(atop("Epiphyte load",paste("(mg cm"^"-2"~")"))))+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        legend.position = "none")
p_epi

# sheath length plot
p_sh <- ggplot()+
  geom_pointrange(data= combo,aes(x=as.factor(order),
                                    y=SheathLengthMean/10,
                                    ymax=(SheathLengthMean+SheathLengthSe)/10,
                                    ymin=(SheathLengthMean-SheathLengthSe)/10,
                                    color=Type,shape=Type))+
  #scale_y_continuous(expand = c(0,10))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","green4"))+
  scale_shape_manual(values=c(8,16))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab("Sheath length (cm)")+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        legend.position="none")
p_sh

# canopy height plot
p_ch <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),
                                    y=CanopyHeightMean/10,
                                    ymax=(CanopyHeightMean+CanopyHeightSe)/10,
                                    ymin=(CanopyHeightMean-CanopyHeightSe)/10,
                                    color=Type,shape=Type))+  
 facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","green4"))+
  scale_shape_manual(values=c(8,16))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab("Canopy height (cm)")+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        legend.position = "none")
p_ch

# blade area plot
p_ba <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),
                                    y=BladeAreaMean,
                                    ymax=BladeAreaMean+BladeAreaSe,
                                    ymin=BladeAreaMean-BladeAreaSe,color=Type,shape=Type))+
  scale_y_continuous(expand = c(0,10))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","green4"))+
  scale_shape_manual(values=c(8,16))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab(expression("Blade area (cm"^"2"~")"))+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        plot.margin = margin(t=3,b=3,r=3,l=5,unit="mm"),
        legend.position = "none")
p_ba

# disease prevalence plot
p_prev <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),
                                    y=PrevalenceMean,
                                    ymax=PrevalenceMean+PrevalenceSe,
                                    ymin=PrevalenceMean-PrevalenceSe,color=Type,shape=Type))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","wheat4"))+
  scale_shape_manual(values=c(8,16))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("Wasting disease prevalence\n(% individuals infected)")+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        plot.margin = margin(t=3,b=3,r=3,l=5,unit="mm"),
        legend.position = "none")
p_prev

# disease severity plot
p_sev <- ggplot()+
  geom_pointrange(data=combo,aes(x=as.factor(order),
                                    y=SeverityMean,
                                    ymax=SeverityMean+SeveritySe,
                                    ymin=SeverityMean-SeveritySe,color=Type,shape=Type))+
  facet_wrap(~Region,nrow = 1,scales = "free_x",strip.position = "bottom")+
  scale_color_manual(values=c("black","wheat4"))+
  scale_shape_manual(values=c(8,16))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  scale_x_discrete(breaks = combo$order,
                   labels = combo$SiteCode)+
  ylab("Wasting disease severity\n(% blade area damaged)")+
  xlab("Site within Region")+
  theme_bw()+
  #  theme(axis.text.x=element_text(angle=90))+
  theme(panel.spacing = unit(0, "lines"),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.border = element_rect(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        strip.text = element_text(size=9),
        plot.margin = margin(t=3,b=3,r=3,l=5,unit="mm"),
        legend.position = "none")
p_sev

# Combine plots into figures ####

# canopy height, sheath length, and epiphyte load for supplemental Figure S5
p_ch1 <- p_ch+xlab(NULL)
p_sh1 <- p_sh+xlab(NULL)

p_all <- plot_grid(
  p_ch1,
  p_sh1,
  p_epi,
  align="v",
  labels=c("A","B","C"),
  ncol=1
)
p_all
ggsave(p_all,filename = "Figures/FigS5_seagrass_metrics.jpg",height=6.75,width=6)
# create high resolution version for paper submission (high-res files not uploaded to github)
ggsave(p_all,filename = "Figures/HighRes/FigS5_seagrass_metrics.tiff",height=6.75,width=6)

# blade area, shoot density, disease prevalence, and disease severity for manuscript Figure 6
p_ba1 <- p_ba+xlab(NULL)
p_den1 <- p_den+xlab(NULL)
p_prev1 <- p_prev+xlab(NULL)

p_fig <- plot_grid(
  p_ba1,
  p_den1,
  p_prev1,
  p_sev,
  align="v",axis="l",
  labels=c("A","B","C","D"),hjust = 0,
  #nrow=3,
  ncol=1
)
p_fig
ggsave(p_fig,filename = "Figures/Fig6_disease_metrics.jpg",height=9,width=6)
ggsave(p_fig,filename = "Figures/HighRes/Fig6_disease_metrics.tiff",height=9,width=6)
