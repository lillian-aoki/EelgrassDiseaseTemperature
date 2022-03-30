# Code files for EelgrassDiseaseTemperature manuscript
# 04_meadow_disease_model

# Last updated 2021-05-20 by Lillian Aoki

# This script uses eelgrass wasting disease survey data and remotely sensed SST data to model effects of temperature anomalies
# and plant and meadow characteristsics on wasting disease prevalence and severity. 

# outputs include Fig S1A and S1B (effect sizes of the prevalence and severity models) and Fig 5 (model predictors + data)

library(tidyverse)
library(betareg)
library(performance)
library(patchwork)
library(cowplot)
library(sjPlot)
region_order <- c("AK","BC","WA","OR","BB","SD")

# Read in and arrange data ####
dat <- read.csv("Data/june19_9y_SST_anomaly_disease.csv")
dat$Region <- ordered(dat$Region,levels=region_order)
dat$CPTempAnomaly <- dat$CDiffMeanHeat
ba <- read.csv("Data/all_survey_metrics_site.csv")
ba$Meadow <- paste(ba$Region,ba$SiteCode,sep="_")
ba <- select(ba,c("Meadow","BladeAreaMean", "DensityShootsMean", "EpiphytePerAreaMean"))
dat <- left_join(dat,ba,by="Meadow")

# Prevalence model ####

# Compare four models, one with cumulative anomaly and blade area as fixed effects, one with cumulative anomaly and shoot density,
# third with only cumulative anomaly and last with all three fixed effects
fit_prev1 <- betareg(PrevalenceMean~CPTempAnomaly + BladeAreaMean,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev2 <- betareg(PrevalenceMean~CPTempAnomaly + DensityShootsMean,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev3 <- betareg(PrevalenceMean~CPTempAnomaly + EpiphytePerAreaMean,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev4 <- betareg(PrevalenceMean~CPTempAnomaly,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev5 <- betareg(PrevalenceMean~ BladeAreaMean,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev6 <- betareg(PrevalenceMean~ DensityShootsMean,
                     data=dat,
                     weights = Count,
                     link = "logit")
fit_prev7 <- betareg(PrevalenceMean~ EpiphytePerAreaMean,
                     data=dat,
                     weights = Count,
                     link = "logit")

df.AIC <- AIC(fit_prev1,fit_prev2,fit_prev3,fit_prev4,fit_prev5,fit_prev6, fit_prev7)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# First model (cumulative anomaly and blade area) is substantially better by AIC.

summary(fit_prev1)
E.prev <- resid(fit_prev1,type='pearson')
F.prev <- fitted(fit_prev1)
plot(E.prev~F.prev)
plot(E.prev~dat$BladeAreaMean)
plot(E.prev~dat$CPTempAnomaly)
# No major issues with the residuals when plotted against covariates and fitted values. Model is acceptable.
# 
# Both cumulative anomaly and blade area are significant. 
# 
## Visualize the best-fitting model effect sizes
prev_names <- c("Leaf\n area", "Cumulative\n SST \nanomaly")
meadow_plot <- plot_model(fit_prev1,
                          type="std",
                          axis.labels = prev_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.2,
                          value.size = 3,
                          axis.lim = c(0.2,15),
                          group.terms = c(1,1))
Sa <- meadow_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.5,2.5))+
  scale_color_manual(values=c("black"))+
  ylab("Scaled estimates of \ndisease prevalence\n odds ratio")+
  labs(tag="A")+
  #xlab("Scaled parameters")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=9))
Sa
# ggsave(filename = "Figures/FigS1A_odds_ratio_meadow_prevalence.jpg", width=4, height=4)
# create high resolution version for submission
# ggsave(filename = "Figures/HighRes/FigS1A_odds_ratio_meadow_prevalence.tiff", width=4, height=4)
# Odds ratio figure is Supplemental Figure S1A (effect sizes of meadow prevalence model)
# Odds ratio plot shows the standardized effect sizes. For a 1 SD increase in Cumulative SST anomaly, the chances of a meadow being completely diseased (prevalence = 100%) effectively doubles. Effect of leaf area is much weaker (only increases prevalence by a factor of 1.14x for each increase in SD).
# 
# Note, SD of Cumulative SST anomaly is 7.4ºC, so for every ~7 degrees of anomaly in the month of June, get a doubling of prevalence
# 
# Visualize model by simulating fit across new data holding one variable constant at a time.

# First simulate change in CPTA
jvaluesPa <- with(dat, seq(from = min(CPTempAnomaly), to = max(CPTempAnomaly), length.out = 100))
# create new data and hold blade area at the median value of the dataset (29.1 cm2)
b.dataPa <- data.frame(CPTempAnomaly=jvaluesPa,BladeAreaMean=median(dat$BladeAreaMean))
predPa <- cbind(
  response=predict(fit_prev1,newdata=b.dataPa,type='response'),
  variance=predict(fit_prev1,newdata=b.dataPa,type='variance'),
  predict(fit_prev1,newdata=b.dataPa,type='quantile',at=c(0.025,0.975)))
preva <- as.data.frame(predPa)
preva <- cbind(preva,b.dataPa)
a <- ggplot(preva,aes(x=CPTempAnomaly))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=PrevalenceMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab("Cumulative positive temperature anomaly (ºC)")+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw()+
  theme(panel.grid = element_blank())

# Next, repeat simulation but hold CPTA constant at median value (11.7 ºC)
jvaluesPb <- with(dat, seq(from = min(BladeAreaMean), to = max(BladeAreaMean), length.out = 100))
# create new data and hold blade area at the median value of the dataset (29.1 cm2)
b.dataPb <- data.frame(BladeAreaMean=jvaluesPb,CPTempAnomaly=median(dat$CPTempAnomaly))
predPb <- cbind(
  response=predict(fit_prev1,newdata=b.dataPb,type='response'),
  variance=predict(fit_prev1,newdata=b.dataPb,type='variance'),
  predict(fit_prev1,newdata=b.dataPb,type='quantile',at=c(0.025,0.975)))
prevb <- as.data.frame(predPb)
prevb <- cbind(prevb,b.dataPb)
b <- ggplot(prevb,aes(x=BladeAreaMean))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=BladeAreaMean,y=PrevalenceMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab(expression("Blade Area (cm"^"2"~")"))+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw()+
  theme(panel.grid = element_blank())
a/b + plot_layout(guides="collect")

# Finally, calculate standardized coefficients (can do this by scaling the inputs to the model?)
SDy <- sd(dat$PrevalenceMean)
SDx1 <- sd(dat$CPTempAnomaly)
SDx2 <- sd(dat$BladeAreaMean)
coefficients(fit_prev1)[2]*SDx1
coefficients(fit_prev1)[3]*SDx2

# Standardized coefficient for CPTA is about 5x greater than the standardized coefficient for blade area.

# Severity model ####

# Repeat same analysis for severity

fit_sev1 <- betareg(SeverityMean~CPTempAnomaly + BladeAreaMean,
                    data=dat,
                    weights = Count,
                    link = "logit")

fit_sev2 <- betareg(SeverityMean~CPTempAnomaly + DensityShootsMean,
                    data=dat,
                    weights = Count,
                    link = "logit")
fit_sev3 <- betareg(SeverityMean~CPTempAnomaly + EpiphytePerAreaMean,
                    data=dat,
                    weights = Count,
                    link = "logit")
fit_sev4 <- betareg(SeverityMean~CPTempAnomaly,
                    data=dat,
                    weights = Count,
                    link = "logit")
fit_sev5 <- betareg(SeverityMean~ BladeAreaMean,
                    data=dat,
                    weights = Count,
                    link = "logit")
fit_sev6 <- betareg(SeverityMean~ DensityShootsMean,
                    data=dat,
                    weights = Count,
                    link = "logit")
fit_sev7 <- betareg(SeverityMean~ EpiphytePerAreaMean,
                    data=dat,
                    weights = Count,
                    link = "logit")
df.AIC <- AIC(fit_sev1,fit_sev2,fit_sev3,fit_sev4,fit_sev5, fit_sev6, fit_sev7)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC

# Note, the model with three fixed predictors is better by AIC for severity. 
# BUT this model is overfitted - the p-values are suddenly extremely small, including for CTempAnomaly. Which doesn't make sense. 
# So we cannot have three predictors on this small of a dataset. 
# 
# Second note, the models with shoot density are worse by AIC than models with blade area. 
# 
# Use the same model as for prevalence (want to know if the same factors drive both).

summary(fit_sev1)
E.sev <- resid(fit_sev1,type='pearson')
F.sev <- fitted(fit_sev1)
plot(E.sev~F.sev)
plot(E.sev~dat$BladeAreaMean)
plot(E.sev~dat$CPTempAnomaly)
# Beach Haven residual is a bit of an outlier perhaps 
# but otherwise no major issues with the residuals when plotted against covariates and fitted values. Model is acceptable.
# 
# Blade area is significant, CPTA is not. 

## Visualize the best-fitting model effect sizes
sev_names <- c("Leaf\n area", "Cumulative\n SST\n anomaly")
meadow_plot2 <- plot_model(fit_sev1,
                           type="std",
                           axis.labels = sev_names,
                           title="",
                           show.p = TRUE,
                           show.values = TRUE,
                           value.size = 3,
                           value.offset = 0.2,
                           axis.lim = c(0.2,15),
                           group.terms = c(1,2))
Sb <- meadow_plot2+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.5,2.5))+
  scale_color_manual(values=c("grey50", "black"))+
  ylab("Scaled estimates of \ndisease severity ratio")+
  #xlab("Scaled parameters")+
  labs(tag="B")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=9))
Sb
# ggsave(filename = "Figures/FigS1B_effect_size_meadow_severity.jpg", width=2.2, height=4)
# create high resolution version for submission (not uploaded)
# ggsave(filename = "Figures/HighRes/FigS1B_effect_size_meadow_severity.tiff", width=4, height=4)

# Visualize model by simulating fit across new data holding one variable constant at a time.

jvaluesSc <- with(dat, seq(from = min(BladeAreaMean), to = max(BladeAreaMean), length.out = 100))
b.dataSc <- data.frame(BladeAreaMean=jvaluesSc,CPTempAnomaly=median(dat$CPTempAnomaly),Region="BC")
predSc <- cbind(
  response=predict(fit_sev1,newdata=b.dataSc,type='response'),
  variance=predict(fit_sev1,newdata=b.dataSc,type='variance'),
  predict(fit_sev1,newdata=b.dataSc,type='quantile',at=c(0.025,0.975)))
sevc <- as.data.frame(predSc)
sevc <- cbind(sevc,b.dataSc)
c <- ggplot(sevc,aes(x=BladeAreaMean))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=BladeAreaMean,y=SeverityMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab(expression("Blade Area (cm"^"2"~")"))+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw()+
  theme(panel.grid = element_blank())

# Don't simulate with CPTA because it's not significant

d <- ggplot()+
  #geom_line(aes(y=response))+
  #geom_line(aes(y=q_0.025),linetype="dashed")+
  #geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=SeverityMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab("Cumulative positive temperature anomaly (ºC)")+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw()+
  theme(panel.grid = element_blank())
(d / c) + plot_layout(guides="collect")

SDy <- sd(dat$PrevalenceMean)
SDx1 <- sd(dat$CPTempAnomaly)
SDx2 <- sd(dat$BladeAreaMean)
coefficients(fit_sev1)[2]*SDx1
coefficients(fit_sev1)[3]*SDx2

# Here, the standardized coeff for blade area is 30x greater in magnitude than the CPTA (which is non-sig anyway).
# 
# Lesion Area model ####
# repeat model selection for lesion area as a disease metric - same predictors as for prevalence and severity
# model is a glm with gamma distribution now (because lesion area is not bounded by 0-1 but rather >0)
fit_les1 <- glm(LesionAreaMean~CPTempAnomaly + BladeAreaMean,
                    data=dat,
                    weights = Count,
                    family = Gamma(link="log"))
fit_les2 <- glm(LesionAreaMean~CPTempAnomaly + DensityShootsMean,
                    data=dat,
                    weights = Count,
                    family = Gamma(link="log"))
fit_les3 <- glm(LesionAreaMean~CPTempAnomaly + EpiphytePerAreaMean,
                data=dat,
                weights = Count,
                family = Gamma(link="log"))
fit_les4 <- glm(LesionAreaMean~CPTempAnomaly,
                    data=dat,
                    weights = Count,
                    family = Gamma(link="log"))
fit_les5 <- glm(LesionAreaMean~ BladeAreaMean,
                    data=dat,
                    weights = Count,
                    family = Gamma(link="log"))
fit_les6 <- glm(LesionAreaMean~ DensityShootsMean,
                    data=dat,
                    weights = Count,
                    family = Gamma(link="log"))
fit_les7 <- glm(LesionAreaMean~ EpiphytePerAreaMean,
                data=dat,
                weights = Count,
                family = Gamma(link="log"))
df.AIC <- AIC(fit_les1,fit_les2,fit_les3,fit_les4,fit_les5,fit_les6, fit_les7)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# use fit_les1, by far the most information by AIC
summary(fit_les1)
E.prev <- resid(fit_les1,type='pearson')
F.prev <- fitted(fit_les1)
plot(E.prev~F.prev)
plot(E.prev~dat$BladeAreaMean)
plot(E.prev~dat$CPTempAnomaly)
plot(E.prev~dat$DensityShootsMean)
# resids are ok

les_names <- c("Leaf\n area", "Cumulative\n SST\n anomaly")
meadow_plot <- plot_model(fit_les3,
                          type="std",
                          axis.labels = les_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.2,
                          value.size = 3,
                          axis.lim = c(0.2,15),
                          group.terms = c(1,2)
)
Sc <- meadow_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.5,2.6))+
  scale_color_manual(values=c("black","grey50"))+
  ylab("Scaled estimates of \nlesion area")+
  labs(tag="C")+
  #xlab("Scaled parameters")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=9))
Sc
# ggsave(filename = "Figures/FigS1C_effect_size_meadow_lesion_area.jpg", width=2, height=4)
# create high resolution version for submission (not uploaded)
# ggsave(filename = "Figures/HighRes/FigS1B_effect_size_meadow_severity.tiff", width=4, height=4)

# Visualize model by simulating fit across new data holding one variable constant at a time.

# First simulate change in CPTA
jvaluesLe <- with(dat, seq(from = min(CPTempAnomaly), to = max(CPTempAnomaly), length.out = 100))
# create new data and hold blade area at the median value of the dataset (29.1 cm2)
b.dataLe <- data.frame(CPTempAnomaly=jvaluesLe,BladeAreaMean=median(dat$BladeAreaMean))
# make predictions on scale of linear predictor
predLe_lin <- predict(fit_les1, newdata=b.dataLe, type="link", se.fit=TRUE)
# use critical interval of 1.96 to calculate 95% for SE
critval <- 1.96
# need to use inverse link function to scale to response variable
predLe <- cbind(
  response=fit_les1$family$linkinv(predLe_lin$fit),
  upper=fit_les1$family$linkinv(predLe_lin$fit + (critval * predLe_lin$se.fit)),
  lower=fit_les1$family$linkinv(predLe_lin$fit - (critval * predLe_lin$se.fit))
)
# predLe <- cbind(
#   response=predict(fit_les1,newdata=b.dataLe,type='response'),
#   #variance=predict(fit_les1,newdata=b.dataLe,type='variance'),
#   se_link=predict(fit_les1,newdata=b.dataLe,type='link',se.fit=TRUE)[[2]],
#   upper=response+in)
lese <- as.data.frame(predLe)
lese <- cbind(lese,b.dataLe)
e <- ggplot(lese,aes(x=CPTempAnomaly))+
  geom_line(aes(y=response))+
  geom_line(aes(y=upper),linetype="dashed")+
  geom_line(aes(y=lower),linetype="dashed")+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=LesionAreaMean,color=Region),size=2)+
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab("Cumulative positive temperature anomaly (ºC)")+
  ylab(expression(paste("Wasting disease \nlesion area (cm"^"2"~")")))+
  theme_bw()+
  theme(panel.grid = element_blank())
e
# Overall, CPTA was significant for determining prevalence but not for severity, which is significantly affected by blade area.

# Plot data and model results for paper ####
a1 <- ggplot(preva,aes(x=CPTempAnomaly))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=PrevalenceMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  guides(color = guide_legend(nrow = 1))+
  xlab("Cumulative positive \ntemperature anomaly (ºC)")+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal"
  )
b1 <-ggplot(prevb,aes(x=BladeAreaMean))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=BladeAreaMean,y=PrevalenceMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab(expression(paste("Leaf area (cm"^2,")")))+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.title = element_blank(),
        legend.position = "")
c1 <- ggplot()+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=SeverityMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab("Cumulative positive \ntemperature anomaly (ºC)")+
  ylab("Wasting disease severity\n (% leaf area damaged)")+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.title = element_blank(),
        legend.position = "")
d1 <- ggplot(sevc,aes(x=BladeAreaMean))+
  geom_line(aes(y=response))+
  geom_line(aes(y=q_0.025),linetype="dashed")+
  geom_line(aes(y=q_0.975),linetype="dashed")+
  geom_point(data=dat,aes(x=BladeAreaMean,y=SeverityMean,color=Region),size=2)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab(expression(paste("Leaf area (cm"^2,")")))+
  ylab("Wasting disease severity\n (% leaf area damaged)")+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.title = element_blank(),
        legend.position = "")
e1 <- ggplot(lese,aes(x=CPTempAnomaly))+
  geom_line(aes(y=response))+
  geom_line(aes(y=upper),linetype="dashed")+
  geom_line(aes(y=lower),linetype="dashed")+
  geom_point(data=dat,aes(x=CPTempAnomaly,y=LesionAreaMean,color=Region),size=2)+
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab("Cumulative positive \ntemperature anomaly (ºC)")+
  ylab(expression(atop("Wasting disease", paste("lesion area (cm"^2,")"))))+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.title = element_blank(),
        legend.position = "")
f1 <- ggplot()+
  geom_point(data=dat,aes(x=BladeAreaMean,y=LesionAreaMean,color=Region),size=2)+
  # scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_color_viridis_d()+
  xlab(expression(paste("Leaf area (cm"^2,")")))+
  ylab(expression(atop("Wasting disease", paste("lesion area (cm"^2,")"))))+
  theme_bw(base_size = 11)+
  theme(panel.grid = element_blank(),
        legend.margin = margin(t=0,r=0,b=0,l=1,unit="pt"),
        legend.key.size = unit(5,unit="mm"),
        legend.title = element_blank(),
        legend.position = "")
## combine with cowplot
a1a <- a1 +theme(legend.position = "")
legend <- get_legend(a1+theme(legend.box.margin = margin(6,0,0,0),
                              legend.direction = "horizontal"))
pcombo <- cowplot::plot_grid(a1a,b1,nrow=1, labels=c("A","B"))
scombo <- cowplot::plot_grid(c1,d1,nrow=1, labels=c("C","D"))
lcombo <- cowplot::plot_grid(e1, f1, nrow=1, labels=c("E", "F"))
total <- cowplot::plot_grid(pcombo,scombo, lcombo, ncol=1)
total_l <- cowplot::plot_grid(total,legend,nrow=2,rel_heights = c(1,.05))
total_l
# output to Fig 5 in the manuscript
ggsave(filename = "Figures/Fig5_meadow_model.jpg", width = 7, height= 9.5)
# create high resolution version
ggsave(filename = "Figures/HighRes/Fig5_meadow_model.tiff", width = 7, height= 9.5)

supp <- cowplot::plot_grid(Sa, Sb, Sc, ggplot()+theme_blank(), nrow=2, rel_widths = c(1,1))
supp
alt <- cowplot::plot_grid(Sa, Sb, Sc, nrow=1, rel_widths = c(1,1,1))
alt
(Sa + Sb) / (Sc+ggplot()+theme_blank())
ggsave(alt, filename = "Figures/FigS1_effect_size_meadow_model.jpg", width = 8, height = 4)
ggsave(filename = "Figures/HighRes/FigS1_effect_size_meadow_model.tiff", width = 8, height = 4)
