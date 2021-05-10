# Code files for EelgrassDiseaseTemperature manuscript
# 05_transect_disease_model

# Last updated 2021-05-08 by Lillian Aoki

# This script uses eelgrass wasting disease survey data and remotely sensed SST data to model effects of temperature anomalies
# and plant and transect characteristsics on wasting disease prevalence and severity. 

# outputs include Fig S2A, S2B S2C, S2D (effect sizes of the prevalence and severity models),
# Fig S6 and S7 (model predictors + data), and Table S6 (transect model comparisons)

library(tidyverse)
library(TMB)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(performance)
library(sjPlot)
region_order <- c("AK","BC","WA","OR","BB","SD")

# Read in data ####
transect <- read.csv("Data/all_survey_metrics_transect.csv")
cpta <- read.csv("Data/june19_9y_SST_anomaly_disease.csv")
cpta <- select(cpta,c("Region","Site","Meadow","CDiffMeanHeat"))
transect <- left_join(transect,cpta,by=c("Region","SiteCode"="Site"))

transect$Region <- ordered(transect$Region,levels=region_order)
transect$Meadow <- as.factor(transect$Meadow)
transect$sBladeAreaMean <- scale(transect$BladeAreaMean,scale=TRUE,center=TRUE)
transect$sDensityShootsMean <- scale(transect$DensityShootsMean,scale=TRUE,center=TRUE)
transect$sEpiphytePerAreaMean <- scale(transect$EpiphytePerAreaMean,scale=TRUE,center=TRUE)
transect$sCDiffMeanHeat <- scale(transect$CDiffMeanHeat,scale=TRUE,center=TRUE)
dat <- transect[-which(is.na((transect$CDiffMeanHeat))),]
dat$Region <- ordered(dat$Region,levels=region_order)
# center and scale predictors for glmm
dat$sBladeAreaMean <- scale(dat$BladeAreaMean,scale=TRUE,center=TRUE)
dat$sDensityShootsMean <- scale(dat$DensityShootsMean,scale=TRUE,center=TRUE)
dat$sEpiphytePerAreaMean <- scale(dat$EpiphytePerAreaMean,scale=TRUE,center=TRUE)
dat$sCDiffMeanHeat <- scale(dat$CDiffMeanHeat,scale=TRUE,center=TRUE)

# Transect-level prevalence model with CPTA n=162 ####
# predictors: Density, Blade Area, Epiphyte mass per area, CPTA, Tidal height, 
# and interactions between fixed effects and tidal height
# random structure: site within region
# 
fit_prev1s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                        sCDiffMeanHeat:TidalHeight+
                        +Region+(1|Meadow),
                      data=dat,
                      weights=CountBlades,
                      family=binomial)
 
# fixed effects are scaled for easier comparison and to avoid convergence issues
# Note, Region as a random effect (1|Region) results in a singular model (variance of Region is 0)
# This suggests the variation between regions is better explained in the residual than as a random effect (or maybe there are
# not enough regions, although 6 should be sufficient)
# fit_prev1s includes Region as a fixed effect
drop1(fit_prev1s)
# by AIC, dropping Region as a fixed effect improves the model by >2, therefore drop Region
fit_prev2s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                        (1|Meadow),
                      data=dat,
                      weights=CountBlades,
                      family=binomial)
# 
### Model selection via AIC, using a priori model set
# remove blade area : tidal height interaction
fit_prev2.1s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                          sDensityShootsMean:TidalHeight+
                          sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)

# remove shoot density : tidal height interaction
fit_prev2.2s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                          sBladeAreaMean:TidalHeight+
                          sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)
# remove epiphyte load : tidal height interaction
fit_prev2.3s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                          sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                          sCDiffMeanHeat:TidalHeight+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)
# remove cumulative anomaly : tidal height interaction
fit_prev2.4s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                          sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                          sEpiphytePerAreaMean:TidalHeight+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)
# no interactions
fit_prev2.5s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)
# no tidal height (or interactions)
fit_prev2.6s <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+
                          (1|Meadow),
                        data=dat,
                        weights=CountBlades,
                        family=binomial)

df.AIC <- data.frame(AIC(fit_prev2s,fit_prev2.1s,fit_prev2.2s,fit_prev2.3s,fit_prev2.4s,
                         fit_prev2.5s,fit_prev2.6s))
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# Based on AIC and model weights, best model is fit_prev2.1s, second model (no blade area:tidal height interaction)
# output of df.AIC goes into Table S6 of the manuscript

# Now need to validate the model and record output
# Assess the model residuals
E2.sim <- simulateResiduals(fit_prev2.1s)
plot(E2.sim)
plot(E2.sim$scaledResiduals~dat$sBladeAreaMean)
plot(E2.sim$scaledResiduals~dat$sDensityShootsMean)
plot(E2.sim$scaledResiduals~dat$sCDiffMeanHeat)
plot(E2.sim$scaledResiduals~as.factor(dat$TidalHeight))
# Some issues with the uniformity of the simulated residuals, 
# but there's not significant dispersion so it is okay to use the model. 

# Look at the model output, performance metrics
summary(fit_prev2.1s)
performance(fit_prev2.1s)
# output from performance goes into Table S2 of the manuscript

## Visualize the best-fitting model output
p1names <- c("Cumulative SST anomaly * Tidal height (upper)", "Epiphyte load * Tidal height (upper)",
             "Shoot density * Tidal height (upper)",
             "Tidal height (upper)","Cumulative SST anomaly","Epiphyte load","Shoot density","Blade area")
p1 <- plot_model(fit_prev2.1s,
                 axis.labels = p1names,
                 title="",
                 show.p = TRUE,
                 show.values = TRUE,
                 value.offset = 0.3,
                 axis.lim = c(0.2,5),
                 group.terms = c(1,2,1,2,1,2,1,2))
a <- p1+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_manual(values=c("grey50","black"))+
  ylab("Scaled estimates of \n disease prevalence odds ratio")+
  #xlab("Scaled parameters")+
  labs(tag="A")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
a
# ggsave(filename = "Figures/FigS2A_effect_size_transect_prevalence.jpg", width=4, height=4)

# Scaled odds ratio estimates indicate the change in odds of disease prevalence, 
# given a change of 1 SD of the predictor. 
# E.g. for a given transect, an increase of 1 SD of the cumulative SST anomaly 
# means the transect is 2.46 times more likely to be diseased. 
#  
# Transect-level prevalence model, for all transects, no CPTA, n=192 ####
transect$Meadow <- paste(transect$Region,transect$SiteCode,sep="_")
fit_prev5 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)

# summary(fit_prev5)
# no issue with Region as a random effect for the full dataset (more balanced?)

### Model selection by AIC
 # remove blade area : tidal height interaction
 fit_prev5.1 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sDensityShootsMean:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)

 # remove shoot density : tidal height interaction
fit_prev5.2 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)
  # remove epiphyte load : tidal height interaction
fit_prev5.3 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)
 # no interactions
 fit_prev5.4 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)
 # no tidal height (or interactions)
 fit_prev5.5 <- glmmTMB(PrevalenceMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+
                       +(1|Region)+(1|Meadow),
                     data=transect,
                     weights=CountBlades,
                     family=binomial)

 df.AIC <- AIC(fit_prev5,fit_prev5.1,fit_prev5.2,fit_prev5.3,fit_prev5.4,fit_prev5.5)
 df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
 df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
 df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
 df.AIC
# By AIC and weight, model 5.1 is the best-fitting (no blade area interaction with tidal height).
# Output of df.AIC goes into Table S6 of the manuscript
 
# Assess the model residuals
E5.sim <- simulateResiduals(fit_prev5.1)
plot(E5.sim)
plot(E5.sim$scaledResiduals~transect$sBladeAreaMean)
plot(E5.sim$scaledResiduals~transect$sDensityShootsMean)
plot(E5.sim$scaledResiduals~as.factor(transect$TidalHeight))

# Look at the model output, performance metrics
summary(fit_prev5.1)
performance(fit_prev5.1)
## Visualize the best-fitting model output
p2names <- c("Epiphyte load * Tidal height (upper)","Shoot density * Tidal height (upper)",
             "Tidal height (upper)","Epiphyte load","Shoot density","Blade area")
p2 <- plot_model(fit_prev5.1,
           axis.labels = p2names,
           title="",
           show.p = TRUE,
           show.values = TRUE,
           value.offset = 0.3,
           axis.lim = c(0.2,10),
           group.terms = c(1,1,1,1,1,1)
           )
c <- p2+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,4))+
  scale_color_manual(values=c("black"))+
  ylab("Scaled estimates of \ndisease prevalence odds ratio")+
  labs(tag = "C")+
  #xlab("Scaled parameters")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
c
#ggsave(filename = "Figures/FigS2C_effect_size_transect_prevalence.jpg", width=4, height = 4)
# Scaled odds ratio estimates indicate the change in odds of disease prevalence, 
# given a change of 1 SD of the predictor. 
# E.g. for a given transect, an increase of 1 SD of the shoot density means the transect is 1.7 times more likely to be diseased.

# Same parameter estimates for the models with and without CPTA (restricted vs not restricted), 
# but restricted model has a higher marginal R2. 
# Therefore, temperature is perhaps important for explaining the prevalence dynamics.

# Transect level severity with CPTA n = 155 ####
# transect-level severity means are calculated using hurdle approach (only diseased blades are included)
# nine transects had zero diseased blades, so exclude those nine for the base severity dataset
sevT <- subset(transect,SeverityMean>0)
sevTdat <- sevT[-which(is.na((sevT$CDiffMeanHeat))),]

# fit with CPTA to start (data = sevTdat)
fit_sev1 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                      sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                     sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Region)+(1|Meadow),
                    data=sevTdat,
                    family=beta_family(link = "logit"))
summary(fit_sev1)
# again Region variance is 0 because of unevenness, can try fitting Region as a fixed effect but this increases the AIC
 fit_sev2 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                       Region+(1|Meadow),
                     data=sevTdat,
                     family=beta_family(link = "logit"))
summary(fit_sev2)
drop1(fit_sev2)
#drop Region from the model
fit_sev3 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Meadow),
                    data=sevTdat,
                    family=beta_family(link = "logit"))
# summary(fit_sev3)
# Model selection via AIC
 # remove blade area : tidal height interaction
 fit_sev3.1s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Meadow),
                      data=sevTdat,
                    family=beta_family(link = "logit"))

 # remove shoot density : tidal height interaction
fit_sev3.2s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeAreaMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Meadow),
                      data=sevTdat,
                    family=beta_family(link = "logit"))
  # remove epiphyte load : tidal height interaction
fit_sev3.3s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sCDiffMeanHeat:TidalHeight+
                      (1|Meadow),
                      data=sevTdat,
                    family=beta_family(link = "logit"))
 # remove cumulative anomaly : tidal height interaction
 fit_sev3.4s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      (1|Meadow),
                      data=sevTdat,
                    family=beta_family(link = "logit"))
 # no interactions
 fit_sev3.5s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       (1|Meadow),
                     data=sevTdat,
                    family=beta_family(link = "logit"))
 # no tidal height (or interactions)
 fit_sev3.6s <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+sCDiffMeanHeat+
                      (1|Meadow),
                     data=sevTdat,
                    family=beta_family(link = "logit"))
 df.AIC <- AIC(fit_sev3,fit_sev3.1s,fit_sev3.2s,fit_sev3.3s,fit_sev3.4s,fit_sev3.5s,fit_sev3.6s)
 df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
 df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
 df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
 df.AIC
# By AIC and weight, the best model is fit_sev3.6 (no tidal height and no interactions)
# output of df.AIC goes into Table S6 of the manuscript
 
# Assess the model residuals
E3.sim <- simulateResiduals(fit_sev3.6s)
plot(E3.sim)
plot(E3.sim$scaledResiduals~sevTdat$sBladeAreaMean)
plot(E3.sim$scaledResiduals~sevTdat$sDensityShootsMean)
plot(E3.sim$scaledResiduals~sevTdat$sEpiphytePerAreaMean)
plot(E3.sim$scaledResiduals~as.factor(sevTdat$TidalHeight))
plot(E3.sim$scaledResiduals~sevTdat$sCDiffMeanHeat)

# Look at the model output, performance metrics
summary(fit_sev3.6s)
performance(fit_sev3.6s)
## Visualize the best-fitting model output
s_cpta_names <- c("Cumulative SST anomaly","Epiphyte load","Shoot density","Blade area")
s_cpta_plot <- plot_model(fit_sev3.6s,
           axis.labels = s_cpta_names,
           title="",
           show.p = TRUE,
           show.values = TRUE,
           value.offset = 0.3,
            #transform = NULL
           #axis.lim = c(0.1,10),

           group.terms = c(1,2,2,2)
           )
b <- s_cpta_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_manual(values=c("black","grey50"))+
  ylab("Scaled estimates of \ndisease severity odds ratio")+
  labs(tag="B")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
b
# ggsave(filename = "Figures/FigS2B_effect_size_transect_severity.jpg", width=4, height=4)

# In contrast to prevalence models, severity model is beta regression, 
# so the scaled parameter estimates estimate the change in the ratio of proportions, 
# rather than the ratio of probabilities (odds). 
# E.g. for an increase in leaf area of 1 SD, the ratio of diseased to non-diseased tissue will change by a factor of 0.70. 
# Only leaf area is significant; effect sizes are smaller and non-significant for other predictors, including cumulative SST anomaly. 
# 
# Transect level severity without CPTA n = 183 ##### 

# fit models using sevT, all transects, no CPTA
fit_sev1 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      (1|Region)+(1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))
summary(fit_sev1)
# Region is also close to zero, try as fixed effect but it increases AIC
fit_sev2 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      Region+(1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))
drop1(fit_sev2)
# drop Region
# basic model for selection
fit_sev4 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))

 # Model selection
 # remove blade area : tidal height interaction
 fit_sev4.1 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sDensityShootsMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))

 # remove shoot density : tidal height interaction
fit_sev4.2 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+
                      (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))
  # remove epiphyte load : tidal height interaction
fit_sev4.3 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                       sBladeAreaMean:TidalHeight+sDensityShootsMean:TidalHeight+
                       (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))
 # no interactions
 fit_sev4.4 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+TidalHeight+
                      (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))

 # no tidal height (or interactions)
 fit_sev4.5 <- glmmTMB(SeverityMean~sBladeAreaMean+sDensityShootsMean+sEpiphytePerAreaMean+
                      (1|Meadow),
                    data=sevT,
                    family=beta_family(link = "logit"))
 df.AIC <- AIC(fit_sev4,fit_sev4.1,fit_sev4.2,fit_sev4.3,fit_sev4.4,
     fit_sev4.5)
 df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
 df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
 df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
 df.AIC
# By AIC and weight, best model is fit_sev4.1, no blade area: tidal height interaction
# Output of df.AIC goes into Table S6 of the manuscript
# Assess the model residuals
E4.sim <- simulateResiduals(fit_sev4.1)
plot(E4.sim)
plot(E4.sim$scaledResiduals~sevT$sBladeAreaMean)
plot(E4.sim$scaledResiduals~sevT$sDensityShootsMean)
plot(E4.sim$scaledResiduals~sevT$sEpiphytePerAreaMean)
plot(E4.sim$scaledResiduals~as.factor(sevT$TidalHeight))

# Look at the model output, performance metrics
summary(fit_sev4.1)
performance(fit_sev4.1)

## Visualize the best-fitting model output
s_all_names <- c("Epiphyte load * Tidal height (upper)","Shoot density * Tidal height (upper)",
                "Tidal height (upper)","Epiphyte load","Shoot density","Blade area")
s_all_plot <- plot_model(fit_sev4.1,
           axis.labels = s_all_names,
           title="",
           show.p = TRUE,
           show.values = TRUE,
           value.offset = 0.3,
            #transform = NULL
           #axis.lim = c(0.1,10),

           group.terms = c(1,1,2,2,1,2)
           )
d <- s_all_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_manual(values=c("black","grey50"))+
  ylab("Scaled estimates of \ndisease severity odds ratio")+
  labs(tag="D")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
d
#ggsave(filename = "Figures/FigS2D_effect_size_transect_severity.jpg", width=4, height=4)
# Scaled estimates indicate change in the ratio of diseased tissue to non-diseased tissue. 
#For a 1 SD increase in shoot density, the ratio of diseased tissue to non-diseased tissue changes by a factor of 1.84. 
# 
# Same estimate for leaf area but for the unrestricted dataset (all transects), shoot density is significant. 
(a + b) / (c + d) + plot_layout(guides="collect")
ggsave(filename = "Figures/FigS2_effect_size_transect.jpg",width=8, height = 8)

# Figure of transect-level correlations ####

## Prevalence ####
p_ba <- ggplot(data=transect,aes(x=BladeAreaMean,y=PrevalenceMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab(expression("Blade area (cm"^"2"~")"))+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

p_den <- ggplot(data=transect,aes(x=DensityShootsMean,y=PrevalenceMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab(expression("Shoot density (shoots m"^"-2"~")"))+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

p_epi <- ggplot(data=transect,aes(x=EpiphytePerAreaMean,y=PrevalenceMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  scale_x_continuous(limits=c(0,0.016))+
  xlab(expression("Epiphyte load per blade area (g cm"^"-2"~")"))+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

p_cpta <- ggplot(data=transect,aes(x=CDiffMeanHeat,y=PrevalenceMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab("Cumulative positive temperature anomaly (ºC)")+
  ylab("Wasting disease prevalence\n (% individuals infected)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

(p_ba + p_den) / (p_epi + p_cpta) +plot_layout(guides="collect")+plot_annotation(tag_levels = c("A","B","C","D"))
ggsave(filename = "Figures/FigS6_transect_regressions_prevalence.jpg",width=10,height=8)

## Severity ####
dat <- subset(transect,SeverityMean>0)

s_ba <- ggplot(data=dat,aes(x=BladeAreaMean,y=SeverityMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab(expression("Blade area (cm"^"2"~")"))+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

s_den <- ggplot(data=dat,aes(x=DensityShootsMean,y=SeverityMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab(expression("Shoot density (shoots m"^"-2"~")"))+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

s_epi <- ggplot(data=dat,aes(x=EpiphytePerAreaMean,y=SeverityMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  scale_x_continuous(limits=c(0,0.016))+
  xlab(expression("Epiphyte load per blade area (g cm"^"-2"~")"))+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

s_cpta <- ggplot(data=dat,aes(x=CDiffMeanHeat,y=SeverityMean,color=Region,shape=TidalHeight))+
  geom_point(size=3)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ 
  scale_shape_discrete(labels=c("Lower","Upper"))+
  xlab("Cumulative positive temperature anomaly (ºC)")+
  ylab("Wasting disease severity\n (% blade area damaged)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

(s_ba + s_den) / (s_epi + s_cpta) +plot_layout(guides="collect")+plot_annotation(tag_levels = c("A","B","C","D"))
ggsave(filename = "Figures/FigS7_transect_regressions_severity.jpg",width=10,height=8)
ggsave(filename = "Figures/HighRes/FigS7_transect_regressions_severity.tiff", width = 10, height = 8)
