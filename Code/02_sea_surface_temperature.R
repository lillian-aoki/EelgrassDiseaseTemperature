# Code files for EelgrassDiseaseTemperature manuscript
# 02_sea_surface_temperature

# Last updated 2021-05-20 by Lillian Aoki

# This script imports and visualizes sea-surface temperature data for the wasting disease survey sites from 2019,
# including exploratory data analysis to identify relevant temperature metrics correlated with wasting disease
# Outputs include Fig 4, Fig S4 and Table S3 in the manuscript

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

june_stress <- june19%>%
  group_by(Region,Site,Meadow) %>%
  mutate(cDiffMean = cumsum(DiffMean),cDiffMeanHeat=cumsum(DiffMeanHeat),
         cDiff90 =cumsum(DiffT90),cDiff90Heat=cumsum(DiffT90Heat))

region <- june19%>%
  group_by(Region,Date)%>%
  mutate(RegMeanMA=mean(Tmean_ma),Reg90MA=mean(T90_ma))

pb <- as.POSIXct(c("2019-06-01","2019-06-15","2019-06-30"))

SST <- ggplot(june_stress,aes(x=Date))+
  geom_line(aes(y=analysed_sst,color=Site))+
  geom_line(data=region,aes(x=Date,y=RegMeanMA),linetype="dashed")+
  facet_wrap(~Region,ncol=1)+
  ylab("Sea surface temperature (ºC)")+
  scale_color_viridis_d()+
  scale_x_datetime(breaks = pb,date_labels = "%b %d" )+
  guides(colour = guide_legend(nrow = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(margin=margin(1,0,1,0,"pt")),
        panel.grid = element_blank(),
        legend.margin = margin(r=0,l=0,unit="mm"),
        strip.text = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        axis.text = element_text(size=9),
        axis.title = element_text(size=11),
        plot.margin=margin(r=10,unit="pt"),
        legend.position = "bottom")
SST
CPTA <- ggplot(june_stress,aes(x=Date))+
  geom_line(aes(y=cDiffMeanHeat,color=Site))+
  facet_wrap(~Region,ncol=1)+
  ylab("Cumulative positive temperature anomaly (ºC)")+
  scale_color_viridis_d()+
  scale_x_datetime(breaks = pb,date_labels = "%b %d" )+
  #guides(colour = guide_legend(nrow = 1))+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(margin=margin(1,0,1,0,"pt")),
        panel.grid = element_blank(),
        legend.margin = margin(r=0,l=0,unit="mm"),
        strip.text = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9),
        axis.text = element_text(size=9),
        axis.title = element_text(size=11),
        plot.margin = margin(5,5,r=10,5,unit = "pt"),
        legend.position="")
CPTA
p1 <- (SST + CPTA )+ plot_layout(guides = 'collect')
p1

SST1 <- SST +theme(legend.position = "")
legend <- get_legend(SST+theme(legend.box.margin = margin(0,0,0,0),
                               legend.direction = "horizontal"))

total <- plot_grid(SST1,CPTA,ncol=2,labels=c("A","B"),hjust = 0)
total_l <- plot_grid(total,legend,nrow=2,rel_heights = c(1,.1))
total_l
# plot of SST and CPTA in June is Fig 4 in the manuscript
ggsave(total_l,filename = "Figures/Fig4_SST_CPTA.jpg",width=4.75,height=6.25)
# create high resolution version for manuscript submission (not uploaded)
ggsave(total_l,filename = "Figures/HighRes/Fig4_SST_CPTA.tiff",width=4.75,height=6.25)

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
ggsave(filename = "Figures/FigS4_daily_SST.jpg",width=6,height=4)
# create high resolution version, not uploaded
ggsave(filename = "Figures/HighRes/FigS4_daily_SST.tiff",width=6,height=4)
