# Code files for EelgrassDiseaseTemperature manuscript
# 02_sea_surface_temperature

# Last updated 2022-04-20 by Lillian Aoki

# This script imports and visualizes sea-surface temperature data for the wasting disease survey sites from 2019,
# including exploratory data analysis to identify relevant temperature metrics correlated with wasting disease
# Outputs include Fig 4, Fig S2, Fig S3 and Table S2 in the manuscript

library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)
library(lubridate)
library(RcppRoll)
region_order <- c("AK","BC","WA","OR","BB","SD")
season_order <- c("Winter","Spring","Summer")

### Read in data ####
# daily SST data from two GHRSST products (MUR and G1SST)
mur <- read.csv("Data/MUR_allsites_9y.csv")
ghr <- read.csv("Data/GHRSST_allsites_combo.csv")
ghr <- dplyr::select(ghr,-SiteCode)
names(ghr)[names(ghr)=="SST"] <- "analysed_sst"
#Combine into a data set with all SST data from 27 sites total
combo <- rbind(mur,ghr)
length(unique(combo$Meadow))
# Exclude MUR temps from before 2011-01-01 (only 9 years of data are available for the G1SST sites)
combo$time <- as.POSIXct(combo$time)
short <- subset(combo,time>"2010-12-31"&time<"2020-01-01")
annual <- subset(short,time>"2018-12-31"&time<"2020-01-01")
annual$Date <- as.Date(annual$Date)
monthly <- annual %>%
  group_by(Region,Site,Meadow,Month=floor_date(Date,unit="month"))%>%
  summarise(meanT=mean(analysed_sst),maxT=max(analysed_sst),rangeT=maxT-min(analysed_sst))
# read in disease data at site level
# disease <- read.csv("Data/all_survey_metrics_site.csv")
# disease$Region <- ordered(disease$Region,levels=region_order)

## input metadata with sample date and disease data
# note this is the same dataset as for disease
meta <- read.csv("Data/all_survey_metrics_site.csv")
meta$Region <- ordered(meta$Region, levels=region_order)
meta <- dplyr::select(meta, c(Year, Region, SiteCode, PrevalenceMean, SeverityMean, LesionAreaMean, SampleDate, SampleJulian))

# input logger data for comparison with SST 
JJA19 <- subset(combo,time>"2019-05-31"&time<"2019-09-01")
hobo <- read.csv("Data//HOBO_JJA_2019.csv")
hobo$DateObs <- as.POSIXct(hobo$DateObs)
hobo_summ <- hobo %>%
  group_by(Region,SiteCode,Day=floor_date(DateObs,unit="day"))%>%
  summarise(TempIS=mean(TempC))

### Raw temperature metrics ####
month_d <- left_join(monthly,meta,by=c("Region","Site"="SiteCode"))
months <- unique(month_d$Month)
rm(corr_table)
corr_table <- list()
# look at correlations between disease levels (prevalence and severity) and temperature metrics
# monthly temperature metrics
for(i in seq_along(months)){
  month=subset(month_d,Month==months[i])
  prev_mean <- cor.test(month$PrevalenceMean,month$meanT)
  prev_max <- cor.test(month$PrevalenceMean,month$maxT)
  prev_range <- cor.test(month$PrevalenceMean,month$rangeT)
  sev_mean <- cor.test(month$SeverityMean,month$meanT)
  sev_max <- cor.test(month$SeverityMean,month$maxT)
  sev_range <- cor.test(month$SeverityMean,month$rangeT)
  les_mean <- cor.test(month$LesionAreaMean, month$meanT)
  les_max <- cor.test(month$LesionAreaMean, month$maxT)
  les_range <- cor.test(month$LesionAreaMean, month$rangeT)
  month_row <- data.frame(Month=months[i],Prev_mean_corr=prev_mean$estimate,Prev_mean_p=prev_mean$p.value,
                          Prev_max_corr=prev_max$estimate,Prev_max_p=prev_max$p.value,
                          Prev_range_corr=prev_range$estimate,Prev_range_p=prev_range$p.value,
                          Sev_mean_corr=sev_mean$estimate,Sev_mean_p=sev_mean$p.value,
                          Sev_max_corr=sev_max$estimate,Sev_max_p=sev_max$p.value,
                          Sev_range_corr=sev_range$estimate,Sev_range_p=sev_range$p.value,
                          Les_mean_corr=les_mean$estimate, Les_mean_p=les_mean$p.value,
                          Les_max_corr=les_max$estimate, Les_max_p=les_max$p.value,
                          Les_range_corr=les_range$estimate, Les_range_p=les_range$p.value)
  corr_table[[i]] <- month_row
}
# bind rows into table that shows Pearson's correlation coefficient and p-value for different temperature metrics
tab1 <- bind_rows(corr_table)
# calculate seasonal temperature metrics
annual$Season[annual$time<"2019-04-01"] <- "Winter"
annual$Season[annual$time>"2019-03-31" & annual$time<"2019-07-01"] <- "Spring"
annual$Season[annual$time>"2019-06-30"&annual$time<"2019-09-01"] <- "Summer"
seasonal <- annual %>%
  group_by(Region,Site,Meadow,Season)%>%
  summarise(meanT=mean(analysed_sst),maxT=max(analysed_sst),rangeT=maxT-min(analysed_sst))

season_d <- left_join(seasonal,meta,by=c("Region","Site"="SiteCode"))
seasons <- na.omit(unique(season_d$Season))
rm(corr_table)
corr_table <- list()
# look at correlation between seasonal temperature metrics and disease
for(i in seq_along(seasons)){
  season=subset(season_d,Season==seasons[i])
  prev_mean <- cor.test(season$PrevalenceMean,season$meanT)
  prev_max <- cor.test(season$PrevalenceMean,season$maxT)
  prev_range <- cor.test(season$PrevalenceMean,season$rangeT)
  sev_mean <- cor.test(season$SeverityMean,season$meanT)
  sev_max <- cor.test(season$SeverityMean,season$maxT)
  sev_range <- cor.test(season$SeverityMean,season$rangeT)
  les_mean <- cor.test(season$LesionAreaMean, season$meanT)
  les_max <- cor.test(season$LesionAreaMean, season$maxT)
  les_range <- cor.test(season$LesionAreaMean, season$rangeT)
  season_row <- data.frame(Season=seasons[i],Prev_mean_corr=prev_mean$estimate,Prev_mean_p=prev_mean$p.value,
                           Prev_max_corr=prev_max$estimate,Prev_max_p=prev_max$p.value,
                           Prev_range_corr=prev_range$estimate,Prev_range_p=prev_range$p.value,
                           Sev_mean_corr=sev_mean$estimate,Sev_mean_p=sev_mean$p.value,
                           Sev_max_corr=sev_max$estimate,Sev_max_p=sev_max$p.value,
                           Sev_range_corr=sev_range$estimate,Sev_range_p=sev_range$p.value,
                           Les_mean_corr=les_mean$estimate, Les_mean_p=les_mean$p.value,
                           Les_max_corr=les_max$estimate, Les_max_p=les_max$p.value,
                           Les_range_corr=les_range$estimate, Les_range_p=les_range$p.value)
  corr_table[[i]] <- season_row
}
tab2 <- bind_rows(corr_table)
# values from tab1 and tab2 go into Table S3
# From tab1 and tab2, no strong correlations between raw temperatures and disease metrics

# Cumulative temperature anomalies ####

# Calculate mean daily temperatures for the combined dataset, using 11-day rolling average
# 11-day window is based on the protocols for identifying marine heat waves (Hobday et al. 2016)
short_summ <- short %>%
  group_by(Meadow,Region,Site,Julian,Lat=lat,Long=lon)%>%
  summarise(Tmean=mean(analysed_sst,na.rm=TRUE),count=length(analysed_sst),T90=quantile(analysed_sst,0.90,na.rm=TRUE))
f11 <- rep(1/11, 11)
short_summ$Tmean_ma <- as.numeric(stats::filter(short_summ$Tmean, f11, sides=2))
short_summ$T90_ma <- as.numeric(stats::filter(short_summ$T90, f11, sides=2))
short_19 <- subset(short,time>"2018-12-31"&time<"2020-01-01")
short_19 <- left_join(short_19,short_summ,by=c("Julian","lat"="Lat","lon"="Long","Meadow","Region","Site"))
short_19$DiffMean <- short_19$analysed_sst-short_19$Tmean_ma
short_19$DiffT90 <- short_19$analysed_sst-short_19$T90_ma

short_19$Region <- ordered(short_19$Region,levels=region_order)
short_19$Date <- as.POSIXct(short_19$time)
short_19$DiffMeanHeat <- short_19$DiffMean
short_19$DiffMeanHeat[short_19$DiffMeanHeat<0] <- 0
short_19$DiffT90Heat <- short_19$DiffT90
short_19$DiffT90Heat[short_19$DiffT90Heat<0] <- 0


# calculate the cumulative anomaly for 30 days prior to sampling
short_19_30 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 30, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 30, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 30, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 30, align = "right",fill=NA))
stress_30 <- short_19_30[,c("Region","Site","Julian","roll_cDiffMeanHeat","roll_cDiff90Heat")]
dat_30 <- left_join(meta,stress_30,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
# plot the 30-day cumulative anomaly vs disease metrics
ggplot(dat_30,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_30,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
ggplot(dat_30,aes(x=roll_cDiffMeanHeat,y=LesionAreaMean,color=Region))+geom_point()
# assess the correlation (Pearson's) between 30-day anomaly and disease metrics
dat_30 <- na.omit(dat_30)
cor.test(dat_30$PrevalenceMean,dat_30$roll_cDiffMeanHeat)
cor.test(dat_30$SeverityMean,dat_30$roll_cDiffMeanHeat)
cor.test(dat_30$LesionAreaMean,dat_30$roll_cDiffMeanHeat)
cor.test(dat_30$PrevalenceMean,dat_30$roll_cDiff90Heat)
cor.test(dat_30$SeverityMean,dat_30$roll_cDiff90Heat)
cor.test(dat_30$LesionAreaMean,dat_30$roll_cDiff90Heat)
# no strong correlation

# Calculate the cumulative anomaly for 60 days prior to sampling
short_19_60 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 60, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 60, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 60, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 60, align = "right",fill=NA))
stress_60 <- short_19_60[,c("Region","Site","Julian","roll_cDiffMeanHeat","roll_cDiff90Heat")]
dat_60 <- left_join(meta,stress_60,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
# plot the 60-day cumulative anomaly vs Prevalence and Severity 
ggplot(dat_60,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_60,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
# assess correlation (Pearson's) between 60=day cumulative anomaly and disease metrics
dat_60 <- na.omit(dat_60)
cor.test(dat_60$PrevalenceMean,dat_60$roll_cDiffMeanHeat)
cor.test(dat_60$SeverityMean,dat_60$roll_cDiffMeanHeat)
cor.test(dat_60$LesionAreaMean,dat_60$roll_cDiffMeanHeat)
cor.test(dat_60$PrevalenceMean,dat_60$roll_cDiff90Heat)
cor.test(dat_60$SeverityMean,dat_60$roll_cDiff90Heat)
cor.test(dat_60$LesionAreaMean,dat_60$roll_cDiff90Heat)
# no strong signal

## Repeat for 14 days (two weeks) prior to sampling
short_19_14 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 14, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 14, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 14, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 14, align = "right",fill=NA))
stress_14 <- short_19_14[,c("Region","Site","Julian","roll_cDiffMeanHeat","roll_cDiff90Heat")]
dat_14 <- left_join(meta,stress_14,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
ggplot(dat_14,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_14,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
dat_14 <- na.omit(dat_14)
cor.test(dat_14$PrevalenceMean,dat_14$roll_cDiffMeanHeat)
cor.test(dat_14$SeverityMean,dat_14$roll_cDiffMeanHeat)
cor.test(dat_14$LesionAreaMean,dat_14$roll_cDiffMeanHeat)
cor.test(dat_14$PrevalenceMean,dat_14$roll_cDiff90Heat)
cor.test(dat_14$SeverityMean,dat_14$roll_cDiff90Heat)
cor.test(dat_14$LesionAreaMean,dat_14$roll_cDiff90Heat)
# no strong signal

# repeat for 90 days prior to sampling
short_19_90 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 90, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 90, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 90, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 90, align = "right",fill=NA))
stress_90 <- short_19_90[,c("Region","Site","Julian","roll_cDiffMeanHeat","roll_cDiff90Heat")]
dat_90 <- left_join(meta,stress_90,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
ggplot(dat_90,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_90,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
dat_90 <- na.omit(dat_90)
cor.test(dat_90$PrevalenceMean,dat_90$roll_cDiffMeanHeat)
cor.test(dat_90$SeverityMean,dat_90$roll_cDiffMeanHeat)
cor.test(dat_90$LesionAreaMean,dat_90$roll_cDiffMeanHeat)
cor.test(dat_90$PrevalenceMean,dat_90$roll_cDiff90Heat)
cor.test(dat_90$SeverityMean,dat_90$roll_cDiff90Heat)
cor.test(dat_90$LesionAreaMean,dat_90$roll_cDiff90Heat)
# no strong signal

# repeat for 45 days prior to sampling
short_19_45 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 45, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 45, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 45, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 45, align = "right",fill=NA))
stress_45 <- short_19_45[,c("Region","Site","Julian","roll_cDiffMeanHeat","roll_cDiff90Heat")]
dat_45 <- left_join(meta,stress_45,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
ggplot(dat_45,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_45,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
dat_45 <- na.omit(dat_45)
cor.test(dat_45$PrevalenceMean,dat_45$roll_cDiffMeanHeat)
cor.test(dat_45$SeverityMean,dat_45$roll_cDiffMeanHeat)
cor.test(dat_45$LesionAreaMean,dat_45$roll_cDiffMeanHeat)
cor.test(dat_45$PrevalenceMean,dat_45$roll_cDiff90Heat)
cor.test(dat_45$SeverityMean,dat_45$roll_cDiff90Heat)
cor.test(dat_45$LesionAreaMean,dat_45$roll_cDiff90Heat)
# no strong signal

# repeat for 21 days prior to sampling
short_19_21 <- short_19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(roll_cDiffMean = roll_sum(DiffMean, 21, align = "right",fill=NA),
         roll_cDiffMeanHeat=roll_sum(DiffMeanHeat, 21, align = "right",fill=NA),
         roll_cDiff90=roll_sum(DiffT90, 21, align = "right",fill=NA),
         roll_cDiff90Heat=roll_sum(DiffT90Heat, 21, align = "right",fill=NA))
stress_21 <- short_19_21[,c("Region","Site","Julian","roll_cDiffMeanHeat")]
dat_21 <- left_join(meta,stress_21,by=c("Region","SiteCode"="Site","SampleJulian"="Julian"))
ggplot(dat_21,aes(x=roll_cDiffMeanHeat,y=PrevalenceMean,color=Region))+geom_point()
ggplot(dat_21,aes(x=roll_cDiffMeanHeat,y=SeverityMean,color=Region))+geom_point()
cor.test(dat_45$PrevalenceMean,dat_21$roll_cDiffMeanHeat)
cor.test(dat_45$SeverityMean,dat_21$roll_cDiffMeanHeat)
cor.test(dat_45$PrevalenceMean,dat_21$roll_cDiff90Heat)
cor.test(dat_45$SeverityMean,dat_21$roll_cDiff90Heat)
# no strong signal

## calculate monthly cumulative anomaly to get at (non-significant) correlations
monthly_anom <- short_19%>%
  group_by(Region,Site,Meadow,Month=floor_date(Date,unit="month")) %>%
  mutate(cDiffMean = cumsum(DiffMean),cDiffMeanHeat=cumsum(DiffMeanHeat),
         cDiff90 =cumsum(DiffT90),cDiff90Heat=cumsum(DiffT90Heat))
month_stress <- monthly_anom %>%
  group_by(Region,Site,Meadow,Month)%>%
  summarise(CDiffMean=sum(DiffMean),CDiffMeanHeat=sum(DiffMean[DiffMean>0]),
            CDiff90=sum(DiffT90),CDiff90Heat=sum(DiffT90[DiffT90>0]))
month_d <- left_join(month_stress,meta,by=c("Region","Site"="SiteCode"))
months <- unique(month_d$Month)

rm(corr_table,try)
corr_table <- list()
for(i in seq_along(months)){
  month=subset(month_d,Month==months[i])
  month <- na.omit(month)
  prev_mean <- cor.test(month$PrevalenceMean,month$CDiffMeanHeat)
  sev_mean <- cor.test(month$SeverityMean,month$CDiffMeanHeat)
  prev_90 <- cor.test(month$PrevalenceMean,month$CDiff90Heat)
  sev_90 <- cor.test(month$SeverityMean,month$CDiff90Heat)
  les_mean <- cor.test(month$LesionAreaMean, month$CDiffMeanHeat)
  les_90 <- cor.test(month$LesionAreaMean, month$CDiff90Heat)
  month_row <- data.frame(Month=months[i],Prev_mean_corr=prev_mean$estimate,Prev_mean_p=prev_mean$p.value,
                          Sev_mean_corr=sev_mean$estimate,Sev_mean_p=sev_mean$p.value,
                          Prev_90_corr=prev_90$estimate, Prev_90_p=prev_90$p.value,
                          Sev_90_corr=sev_90$estimate, Sev_90_p=sev_90$p.value,
                          Les_mean_corr=les_mean$estimate, Les_mean_p=les_mean$p.value,
                          Les_90_corr=les_90$estimate, Les_90_p=les_90$p.value)
  corr_table[[i]] <- month_row
}
#
try <- bind_rows(corr_table)

# After data exploration, the strongest correlation between temperature metrics and diseae is the cumulative anomaly in June
# Both the anomaly above the mean and the anomaly above the 90th percentile are similar in magnitude and significance
# But, the cumulative anomaly above the mean shows a more consistent relationship across all the sites 
# (versus being driven by a few very hot sites). 
# Therefore, select cumulative anomaly in the month of June as the temperature metric/predictor for further analysis


# Plot June SST and CPTA ####
june19 <- subset(short_19,Julian>151&Julian<182)
june19$RegionName <- plyr::revalue(june19$Region, c("AK"="Alaska", "BC"="British Columbia",
                                              "WA" = "Washington", "OR" = "Oregon", 
                                              "BB"="California - Bodega Bay", "SD" = "California - San Diego"))

june_stress <- june19%>%
  group_by(Region, RegionName,Site,Meadow) %>%
  mutate(cDiffMean = cumsum(DiffMean),cDiffMeanHeat=cumsum(DiffMeanHeat),
         cDiff90 =cumsum(DiffT90),cDiff90Heat=cumsum(DiffT90Heat))

region <- june19%>%
  group_by(Region, RegionName,Date)%>%
  mutate(RegMeanMA=mean(Tmean_ma),Reg90MA=mean(T90_ma))

pb <- as.POSIXct(c("2019-06-01","2019-06-15","2019-06-30"))

SST <- ggplot(june_stress,aes(x=Date))+
  geom_line(aes(y=analysed_sst,color=Site))+
  geom_line(data=region,aes(x=Date,y=RegMeanMA),linetype="dashed")+
  facet_wrap(~RegionName,ncol=1)+
  ylab("Sea surface temperature (ºC)")+
  scale_color_viridis_d()+
  scale_x_datetime(breaks = pb,date_labels = "%b %d" )+
  guides(colour = guide_legend(nrow = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(margin=margin(1,0,1,0,"pt")),
        panel.grid = element_blank(),
        legend.margin = margin(r=0,l=0,unit="mm"),
        strip.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        plot.margin = margin(5,5,r=10,5,unit = "pt"),
        legend.position = "bottom")
SST
CPTA <- ggplot(june_stress,aes(x=Date))+
  geom_line(aes(y=cDiffMeanHeat,color=Site))+
  facet_wrap(~RegionName,ncol=1)+
  ylab("Cumulative positive temperature anomaly (ºC)")+
  scale_color_viridis_d()+
  scale_x_datetime(breaks = pb,date_labels = "%b %d" )+
  #guides(colour = guide_legend(nrow = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(margin=margin(1,0,1,0,"pt")),
        panel.grid = element_blank(),
        legend.margin = margin(r=0,l=0,unit="mm"),
        strip.text = element_text(size=10),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        plot.margin = margin(5,5,r=10,5,unit = "pt"),
        legend.position="")
CPTA
p1 <- (SST + CPTA )+ plot_layout(guides = 'collect')
p1

SST1 <- SST +theme(legend.position = "")
CPTA1 <- CPTA+theme(legend.position = "")
legend <- get_legend(SST+theme(legend.box.margin = margin(0,0,0,0),
                               legend.direction = "horizontal"))

total <- plot_grid(SST1,CPTA1,ncol=2,labels=c("a","b"),hjust = 0)
total_l <- plot_grid(total,legend,nrow=2,rel_heights = c(1,.1))
total_l
# plot of SST and CPTA in June is Fig 4 in the manuscript
ggsave(total_l,filename = "Figures/Fig4_SST_CPTA.jpg",width=6.5,height=6.25)
# create high resolution version for manuscript submission (not uploaded)
ggsave(total_l,filename = "Figures/HighRes/Fig4_SST_CPTA.tiff",width=6.5,height=6.25)

# Plot seasonal SST ranges ####
short$Month <- floor_date(short$time,unit="month")
#short$Month <- as.POSIXct.Date(short$Month,format="%B")
short$Month.w <- as.character.Date(short$Month,format="%B")
short$Season[short$Month.w=="January"|short$Month.w=="February"|short$Month.w=="March"] <- "Winter"
short$Season[short$Month.w=="April"|short$Month.w=="May"|short$Month.w=="June"] <- "Spring"
short$Season[short$Month.w=="July"|short$Month.w=="August"] <- "Summer"

seasonal_region <- short[-which(is.na(short$Season)),] %>%
  group_by(Region,Season,Date)%>%
  summarise(DailyTemp=mean(analysed_sst))
seasonal_region$Region <- ordered(seasonal_region$Region,levels=region_order)
seasonal_region$Season <- ordered(seasonal_region$Season,levels=season_order)
ggplot(seasonal_region,aes(x=Region,y=DailyTemp,fill=Region))+geom_boxplot()+
  facet_wrap(~Season)+
  ylab("Daily SST (ºC)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=11),
        axis.text = element_text(size=10),
        strip.text = element_text(size=10))
ggsave(filename = "Figures/FigS2_daily_SST.jpg",width=6,height=4)
# create high resolution version, not uploaded
ggsave(filename = "Figures/HighRes/FigS2_daily_SST.tiff",width=6,height=4)

# Plot logger temps vs SST temps ####
# Note the comparison is for Aug 2019 because the only summer month with good logger coverage and SST for all sites
RS_IS <- right_join(JJA19,hobo_summ,by=c("Region","Site"="SiteCode","time"="Day"))
noRS <- which(is.na(RS_IS$analysed_sst))
noIS <- which(is.na(RS_IS$TempIS))
RS_IS <- na.omit(RS_IS)
RS_IS <- subset(RS_IS,Meadow!="BC_B")
RS_IS_aug <- subset(RS_IS,time>"2019-07-31"&time<"2019-09-01")
RS_IS_aug$Region <- ordered(RS_IS_aug$Region, levels=region_order)
RS_IS_aug$Meadow <- paste(RS_IS_aug$Region, RS_IS_aug$Site, sep="_")
RS_IS_aug$Meadow <- ordered(RS_IS_aug$Meadow, levels=c("AK_A", "AK_B", "AK_D", "AK_E", "AK_F",
                                                       "BC_A", "BC_C", "BC_D", 
                                                       "WA_A", "WA_B", "WA_C", "WA_D", "WA_E",
                                                       "OR_B", "OR_C", "OR_D", "OR_E",
                                                       "BB_A", "BB_B", "BB_C", "BB_D", "BB_E",
                                                       "SD_A"))
# can create one plot for all data but it's a little messy
all <- ggplot(data=RS_IS_aug,aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  #facet_wrap(~Meadow, scales= "free_x")+
  facet_grid(rows=vars(Region),cols = vars(Site), scales = "free_x")+
  # facet_wrap(Region~Site, scales = "free_x")+
  scale_color_viridis_d(begin = 0, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  xlab("Remotely sensed SST (ºC)")+
  ylab("In situ daily mean temp (ºC)")+
  # labs(title="Comparing satellite and in situ temperatures - OR",
  #      subtitle = "JJA 2019")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.background = element_blank())
# ggsave(filename = "Figures/FigS3_logger_sst_correlations_aug.jpg", width = 6.5, height=7)
# ggsave(filename = "Figures/HighRes/FigS3_logger_sst_correlations_aug.tiff", width = 6.5, height=7)

# instead break out into rows and then recombine with patchwork()
ak <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="AK",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1)+
  scale_color_viridis_d(begin = 0, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none")
ak
ak_row <- plot_grid(ak, NULL, ncol = 2, rel_widths = c(5,0))
ak_row
bc <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="BC",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1, ncol = 5)+
  scale_color_viridis_d(begin = 0.2, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  scale_y_continuous(limits=c(12.5, 16.5), breaks=c(13,14,15,16))+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
bc
bc_row <- plot_grid(bc, NULL, ncol = 2, rel_widths = c(3,2))

wa <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="WA",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1, ncol = 5)+
  scale_color_viridis_d(begin = 0.4, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
wa
wa_row <- plot_grid(wa, NULL, ncol = 2, rel_widths = c(5,0))

or <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="OR",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1, ncol = 5)+
  scale_color_viridis_d(begin = 0.6, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
or
or_row <- plot_grid(or, NULL, ncol = 2, rel_widths = c(4,1))

bb <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="BB",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1, ncol = 5)+
  scale_color_viridis_d(begin = 0.8, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
bb
bb_row <- plot_grid(bb, NULL, ncol = 2, rel_widths = c(5,0))

sd <- ggplot(data=RS_IS_aug[RS_IS_aug$Region=="SD",],aes(x=analysed_sst,y=TempIS, color=Region))+geom_point()+
  stat_smooth(method="lm",col="dark grey")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  facet_wrap(~Meadow, scales= "free_x", nrow = 1, ncol = 5)+
  scale_color_viridis_d(begin = 1, end = 1)+
  scale_x_continuous(breaks = waiver(), n.breaks = 3)+
  scale_y_continuous(limits=c(25, 27), breaks=c(25,26,27))+
  xlab(NULL)+
  ylab(NULL)+
  # xlab("Remotely sensed SST (ºC)")+
  # ylab("In situ daily mean temp (ºC)")+
  theme_bw(base_size = 12)+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
sd
sd_row <- plot_grid(sd, NULL, ncol = 2, rel_widths = c(1,4)) +plot_layout(guides="collect")

legend <- get_legend(all)

y_label <- ggplot()+
  geom_text(aes(x=1, y=1, label="In situ daily mean temperature (ºC)"), angle=90)+
  theme_void()
x_label <- ggplot()+
  geom_text(aes(x=1, y=1, label="Remotely sensed SST (ºC)"))+
  theme_void()+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing= unit(c(0, 0, 0, 0), "cm")
        )

middle <- ((ak_row / bc_row / wa_row / or_row / bb_row / sd_row)/x_label) + 
  plot_layout(nrow = 7, heights = c(1,1,1,1,1,1,0.25))
middle
full <- (y_label | middle | legend)+plot_layout(ncol =3, widths = c(0.2,5,0.5))
full
# (y_label | ((ak_row / bc_row / wa_row / or_row / bb_row / sd_row)/x_label) | legend) + 
#   plot_layout(ncol =3, widths = c(0.2,5,0.5), nrow = 2, heights = c(20, 0.1))
ggsave("Figures/FigS3_logger_sst_correlations_aug.jpg",width = 6.5, height=7)
ggsave("Figures/HighRes/FigS3_logger_sst_correlations_aug.tiff", width = 6.5, height = 7)
