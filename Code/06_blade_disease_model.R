# Code files for EelgrassDiseaseTemperature manuscript
# 06_blade_disease_model

# Last updated 2021-05-20 by Lillian Aoki

# This script uses eelgrass wasting disease survey data and remotely sensed SST data to model effects of temperature anomalies
# and plant and blade characteristsics on wasting disease prevalence and severity. 

# outputs include Fig S6 (effect sizes of the prevalence, lesion area, and severity models) and 
# Table S6-S7 (blade model comparisons)

library(tidyverse)
library(TMB)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(performance)
library(sjPlot)

# Read in data ###

bladeWD <- read.csv("Data/disease_2019.csv")
bladeWD <- bladeWD[,c("Prevalence","Severity","LesionArea","BladeArea","TidalHeight","SampleId","Region","SiteCode","Transect")]
epi <- read.csv("Data/all_survey_metrics_transect.csv")
epi <- epi[,c("Region","SiteCode","TidalHeight","Transect","EpiphytePerAreaMean","DensityShootsMean")]
bladeWD <- left_join(bladeWD,epi,by=c("Region","SiteCode","TidalHeight","Transect"))
cpta <- read.csv("Data/june19_9y_SST_anomaly_disease.csv")
cpta <- select(cpta,c("Region","Site","Meadow","CDiffMeanHeat"))
bladeWD <- left_join(bladeWD,cpta,by=c("Region","SiteCode"="Site"))
bladeWD$sBladeArea <- scale(bladeWD$BladeArea,center=TRUE,scale=TRUE)
bladeWD$sDensityShoots <- scale(bladeWD$DensityShoots,center=TRUE,scale=TRUE)
bladeWD$sEpiphytePerAreaMean <- scale(bladeWD$EpiphytePerAreaMean,center=TRUE,scale=TRUE)
bladeWD$sCDiffMeanHeat <- scale(bladeWD$CDiffMeanHeat,center=TRUE,scale=TRUE)
bladeWD$MeadowId <- paste(bladeWD$Region,bladeWD$SiteCode,sep=".")
bladeWD$TransectId <- paste(bladeWD$MeadowId,bladeWD$Transect,sep="_")
bladeWD$Region <- as.factor(bladeWD$Region)
dat <- bladeWD[-which(is.na((bladeWD$CDiffMeanHeat))),]
#rescale continuous predictors based on the restricted dataset
dat$sBladeArea <- scale(dat$BladeArea,center=TRUE,scale=TRUE)
dat$sDensityShoots <- scale(dat$DensityShoots,center=TRUE,scale=TRUE)
dat$sEpiphytePerAreaMean <- scale(dat$EpiphytePerAreaMean,center=TRUE,scale=TRUE)
dat$sCDiffMeanHeat <- scale(dat$CDiffMeanHeat,center=TRUE,scale=TRUE)

# Blade level prevalence with CPTA (n=3177) ####
fit_prev1 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                    sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                       (1|Region)+(1|MeadowId)+(1|TransectId),
                    data=dat,
                    family=binomial(link="logit"))
summary(fit_prev1)
# so region again isn't the source of variation
# refit with Region as a fixed effect
fit_prev2 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                    sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                       Region+(1|MeadowId)+(1|TransectId),
                     data=dat,
                     family=binomial(link="logit"))
drop1(fit_prev2)
# dropping region improves AIC marginally
fit_prev3 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                       sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                       sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                       (1|MeadowId)+(1|TransectId),
                     data=dat,
                     family=binomial(link="logit"))
 
### Model selection by AIC
# remove blade area : tidal height interaction
fit_prev3.1 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                         sDensityShoots:TidalHeight+
                         sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))

# remove shoot density : tidal height interaction
fit_prev3.2 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                         sBladeArea:TidalHeight+
                         sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))
# remove epiphyte load : tidal height interaction
fit_prev3.3 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                         sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                         sCDiffMeanHeat:TidalHeight+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))
# remove cumulative anomaly : tidal height interaction
fit_prev3.4 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                         sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                         sEpiphytePerAreaMean:TidalHeight+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))
# no interactions
fit_prev3.5 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))

# no tidal height (or interactions)
fit_prev3.6 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+
                         (1|MeadowId)+(1|TransectId),
                       data=dat,
                       family=binomial(link="logit"))
df.AIC <- AIC(fit_prev3,fit_prev3.1,fit_prev3.2,fit_prev3.3,fit_prev3.4,
              fit_prev3.5,fit_prev3.6)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC

# By AIC and weight, fit_prev3.3 (no epiphyte load: tidal height interaction) is best.
# 
# # Assess the model residuals 
E3.sim <- simulateResiduals(fit_prev3.3)
plot(E3.sim)
plot(E3.sim$scaledResiduals~dat$sBladeArea)
plot(E3.sim$scaledResiduals~dat$sDensityShoots)
plot(E3.sim$scaledResiduals~dat$sEpiphytePerAreaMean)
plot(E3.sim$scaledResiduals~as.factor(dat$TidalHeight))
plot(E3.sim$scaledResiduals~dat$sCDiffMeanHeat)

# # Look at the model output, performance metrics
summary(fit_prev3.3)
performance(fit_prev3.3)
## Visualize the best-fitting model output
p_cpta_names <- c("Cumulative SST anomaly * Tidal height (upper)",
                  "Shoot density * Tidal height (upper)", "Leaf area * Tidal height (upper)",
                  "Tidal height (upper)", "Cumulative SST anomaly","Epiphyte load","Shoot density","Leaf area")
p_cpta_plot <- plot_model(fit_prev3.3,
                          axis.labels = p_cpta_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.4,
                          #transform = NULL
                          #axis.lim = c(0.1,10),

                          group.terms = c(1,1,2,1,2,2,1,2))
a <- p_cpta_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \ndisease prevalence odds ratio")+
  labs(tag="A")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
a
# Blade level prevalence, no CPTA n = 3702 #### 
# This is the prevalence model without any temperature terms, with the full 3702 blade dataset.
# For this model, Region can be included in the random part

fit_prev5 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                       sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                       (1|Region)+(1|MeadowId)+(1|TransectId),
                     data=bladeWD,
                     family=binomial(link="logit"))
summary(fit_prev5)


### Model selection by AIC
# remove blade area : tidal height interaction
fit_prev5.1 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                         sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                         (1|Region)+(1|MeadowId)+(1|TransectId),
                       data=bladeWD,
                       family=binomial(link="logit"))

# remove shoot density : tidal height interaction
fit_prev5.2 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                         sBladeArea:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                         (1|Region)+(1|MeadowId)+(1|TransectId),
                       data=bladeWD,
                       family=binomial(link="logit"))
# remove epiphyte load : tidal height interaction
fit_prev5.3 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                         sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                         (1|Region)+(1|MeadowId)+(1|TransectId),
                       data=bladeWD,
                       family=binomial(link="logit"))
# no interactions
fit_prev5.4 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                         (1|Region)+(1|MeadowId)+(1|TransectId),
                       data=bladeWD,
                       family=binomial(link="logit"))
# no tidal height (or interactions)
fit_prev5.5 <- glmmTMB(Prevalence~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+
                         (1|Region)+(1|MeadowId)+(1|TransectId),
                       data=bladeWD,
                       family=binomial(link="logit"))
df.AIC <- AIC(fit_prev5,fit_prev5.1,fit_prev5.2,fit_prev5.3,fit_prev5.4,
              fit_prev5.5)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC

# By AIC and weight, best-fitting model is the full model (fit_prev5).

# Assess the model residuals
E5.sim <- simulateResiduals(fit_prev5)
plot(E5.sim)
plot(E5.sim$scaledResiduals~bladeWD$sBladeArea)
plot(E5.sim$scaledResiduals~bladeWD$sDensityShoots)
plot(E5.sim$scaledResiduals~bladeWD$sEpiphytePerAreaMean)
plot(E5.sim$scaledResiduals~as.factor(bladeWD$TidalHeight))

# Look at the model output, performance metrics
summary(fit_prev5)
performance(fit_prev5)
## Visualize the best-fitting model output
p_all_names <- c("Epiphyte load * Tidal height (upper)", "Shoot density * Tidal height (upper)",
                 "Blade area * Tidal height (upper)","Tidal height (upper)",
                 "Epiphyte load","Shoot density","Leaf area")
p_all_plot <- plot_model(fit_prev5,
                         axis.labels = p_all_names,
                         title="",
                         show.p = TRUE,
                         show.values = TRUE,
                         value.offset = 0.4,
                         #transform = NULL
                         #axis.lim = c(0.1,10),
                         group.terms = c(1,1,2,2,1,1,2))
b <- p_all_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \ndisease prevalence odds ratio")+
  labs(tag="B")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
b

# Blade level severity, with CPTA (n=1573) ####
# Severity at the blade level is the second part of the hurdle, so it is only modeled on the diseased blades.
# Same approach to having two models, one with fewer replicates because of the availability of SST for the anomaly calculation.
# These models are fit with beta regression, as severity is a non-count proportion.

sev <- subset(bladeWD,Prevalence==1)
sev$sBladeArea <- scale(sev$BladeArea,center=TRUE,scale=TRUE)
sev$sDensityShoots <- scale(sev$DensityShoots,center=TRUE,scale=TRUE)
sev$sEpiphytePerAreaMean <- scale(sev$EpiphytePerAreaMean,center=TRUE,scale=TRUE)
sev$sCDiffMeanHeat <- scale(sev$CDiffMeanHeat,center=TRUE,scale=TRUE)
sevdat <- sev[-which(is.na((sev$CDiffMeanHeat))),]
sevdat$sBladeArea <- scale(sevdat$BladeArea,center=TRUE,scale=TRUE)
sevdat$sDensityShoots <- scale(sevdat$DensityShoots,center=TRUE,scale=TRUE)
sevdat$sEpiphytePerAreaMean <- scale(sevdat$EpiphytePerAreaMean,center=TRUE,scale=TRUE)
sevdat$sCDiffMeanHeat <- scale(sevdat$CDiffMeanHeat,center=TRUE,scale=TRUE)

# use beta regression in glmmTMB

fit_sev1 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                    sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Region)+(1|MeadowId)+(1|TransectId),
                               data=sevdat,
                               family=beta_family(link = "logit"))
summary(fit_sev1)
# Region doesn't work as a random effect, try as fixed
fit_sev2 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                    sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      Region+(1|MeadowId)+(1|TransectId),
                    data=sevdat,
                    family=beta_family(link = "logit"))
summary(fit_sev2)
drop1(fit_sev2)
# drop Region as a fixed effect
fit_sev3 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                      sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|MeadowId)+(1|TransectId),
                    data=sevdat,
                    family=beta_family(link = "logit"))

### Model selection by AIC
# remove blade area : tidal height interaction
fit_sev3.1 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sDensityShoots:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))

# remove shoot density : tidal height interaction
fit_sev3.2 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))
# remove epiphyte load : tidal height interaction
fit_sev3.3 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))
# remove cumulative anomaly : tidal height interaction
fit_sev3.4 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))
# no interactions
fit_sev3.5 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))
# no tidal height (or interactions)
fit_sev3.6 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                      family=beta_family(link = "logit"))
df.AIC <- AIC(fit_sev3,fit_sev3.1,fit_sev3.2,fit_sev3.3,fit_sev3.4,
              fit_sev3.5,fit_sev3.6)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC

# By AIC and weight, best model is fit_sev3.4 (no cumulative anomlay:tidal height interaction)

# # Assess the model residuals 
Es.sim <- simulateResiduals(fit_sev3.4)
plot(Es.sim)
plot(Es.sim$scaledResiduals~sevdat$sBladeArea)
plot(Es.sim$scaledResiduals~sevdat$sDensityShoots)
plot(Es.sim$scaledResiduals~sevdat$sEpiphytePerAreaMean)
plot(Es.sim$scaledResiduals~as.factor(sevdat$TidalHeight))
plot(Es.sim$scaledResiduals~sevdat$sCDiffMeanHeat)
testDispersion(Es.sim)
# # Some issues with dispersion here, but we don't need to worry too much for a beta distribution
# Adding a dispersion formula (variable dispersion) doesn't improve the diagnostics
# The residuals plotted against model covariate don't show patterns and have pretty good uniformity, so keep the 3.4 model
# fit_sev_disp <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
#                       sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
#                         sEpiphytePerAreaMean:TidalHeight+
#                       (1|MeadowId)+(1|TransectId),
#                       dispformula = ~sBladeArea,
#                     data=sevdat,
#                     family=beta_family(link = "logit"))
# Edisp.sim <- simulateResiduals(fit_sev_disp)
# plot(Edisp.sim)
# summary(fit_sev_disp)
# 
# # Look at the model output, performance metrics
summary(fit_sev3.4)
performance(fit_sev3.4)
# 
# ## Visualize the best-fitting model output
s_cpta_names <- c("Epiphyte load * Tidal height (upper)",
                  "Shoot density * Tidal height (upper)", "Leaf area * Tidal height (upper)",
                  "Tidal height (upper)", "Cumulative SST anomaly","Epiphyte load","Shoot density","Leaf area")
s_cpta_plot <- plot_model(fit_sev3.4,
                          axis.labels = s_cpta_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.4,
                          #transform = NULL
                          #axis.lim = c(0.1,10),
                          group.terms = c(1,1,1,2,2,2,1,2))
e <- s_cpta_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \ndisease severity odds ratio")+
  labs(tag="E")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
e
# Blade level severity, no CPTA, n=1853 ####
# Follow up model includes all diseased blades, no temperature anomaly term.
# Also fit with beta regression.
# 
fit_sev1 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                      sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                      (1|Region)+(1|MeadowId)+(1|TransectId),
                    data=sev,
                    family=beta_family(link = "logit"))
summary(fit_sev1)
performance(fit_sev1)
drop1(fit_sev1)
# Region does work as a random effect
fit_sev2 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                      sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                      (1|Region)+(1|MeadowId)+(1|TransectId),
                    data=sev,
                    family=beta_family(link = "logit"))
summary(fit_sev2)
performance(fit_sev2)
icc(fit_sev2,by_group = TRUE)


### Model selection by AIC
# remove blade area : tidal height interaction
fit_sev2.1 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=beta_family(link = "logit"))

# remove shoot density : tidal height interaction
fit_sev2.2 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sBladeArea:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=beta_family(link = "logit"))
# remove epiphyte load : tidal height interaction
fit_sev2.3 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=beta_family(link = "logit"))
# no interactions
fit_sev2.4 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=beta_family(link = "logit"))
# no tidal height (or interactions)
fit_sev2.5 <- glmmTMB(Severity~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=beta_family(link = "logit"))
df.AIC <- AIC(fit_sev2,fit_sev2.1,fit_sev2.2,fit_sev2.3,fit_sev2.4,
              fit_sev2.5)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# By AIC and weight, best model is 2.1, no blade area:tidal height interaction
 
# Assess the model residuals
E2.sim <- simulateResiduals(fit_sev2.1)
plot(E2.sim)
plot(E2.sim$scaledResiduals~sev$sBladeArea)
plot(E2.sim$scaledResiduals~sev$sDensityShoots)
plot(E2.sim$scaledResiduals~sev$sEpiphytePerAreaMean)
plot(E2.sim$scaledResiduals~as.factor(sev$TidalHeight))
# Again as with blade severity cpta model, there's an issue with underdispersion, but adding a variable dispersion formula
# doesn't improve residual diagnostics, and the residuals look fine plotted against model covariates. So, accept model.
# Look at the model output, performance metrics
summary(fit_sev2.1)
performance(fit_sev2.1)
## Visualize the best-fitting model output
s_all_names <- c("Epiphyte load * Tidal height (upper)", "Shoot density * Tidal height (upper)",
                 "Tidal height (upper)",
                 "Epiphyte load","Shoot density","Leaf area")
s_all_plot <- plot_model(fit_sev2.1,
                         axis.labels = s_all_names,
                         title="",
                         show.p = TRUE,
                         show.values = TRUE,
                         value.offset = 0.4,
                         #transform = NULL
                         #axis.lim = c(0.1,10),
                         group.terms = c(1,1,1,2,1,1))
f <- s_all_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \ndisease severity odds ratio")+
  labs(tag="F")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
f

# Blade level lesion area, with CPTA n = 1573 ####
# Use the same dataset as for severity.
# Model is a glmm with Gamma distribution
fit_les <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                      sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                      (1|Region)+(1|MeadowId)+(1|TransectId),
                    data=sevdat,
                    family=Gamma(link="log"))
summary(fit_les)
# can include Region as a random effect here
# Model selection by AIC
# remove blade area : tidal height interaction
fit_les1 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sDensityShoots:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))

# remove shoot density : tidal height interaction
fit_les2 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))
# remove epiphyte load : tidal height interaction
fit_les3 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        sCDiffMeanHeat:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))
# remove cumulative anomaly : tidal height interaction
fit_les4 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        sEpiphytePerAreaMean:TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))
# no interactions
fit_les5 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+TidalHeight+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))
# no tidal height (or interactions)
fit_les6 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+sCDiffMeanHeat+
                        (1|MeadowId)+(1|TransectId),
                      data=sevdat,
                    family=Gamma(link="log"))
df.AIC <- AIC(fit_les,fit_les1,fit_les2,fit_les3,fit_les4,
              fit_les5,fit_les6)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# Based on AIC and weight, use fit_les3, no epiphyte load:tidal height interaction
# # Assess the model residuals 
El.sim <- simulateResiduals(fit_les3)
plot(El.sim)
plot(El.sim$scaledResiduals~sevdat$sBladeArea)
plot(El.sim$scaledResiduals~sevdat$sDensityShoots)
plot(El.sim$scaledResiduals~sevdat$sEpiphytePerAreaMean)
plot(El.sim$scaledResiduals~as.factor(sevdat$TidalHeight))
plot(El.sim$scaledResiduals~sevdat$sCDiffMeanHeat)
# Not too worried about non-normality, with thousands of observations. Dispersion is okay, so accept the model
summary(fit_les3)
performance(fit_les3)
# Visualize effect sizes
l_cpta_names <- c("Cumulative anomaly * Tidal height (upper)",
                  "Shoot density * Tidal height (upper)", "Leaf area * Tidal height (upper)",
                  "Tidal height (upper)", "Cumulative SST anomaly","Epiphyte load","Shoot density","Leaf area")
l_cpta_plot <- plot_model(fit_les3,
                          axis.labels =l_cpta_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.4,
                          #transform = NULL
                          #axis.lim = c(0.1,10),
                          group.terms = c(1,1,2,1,2,1,1,2))
c <- l_cpta_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \nlesion area")+
  labs(tag="C")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
c

# Blade level lesion area no CPTA, n = 1853 ####
# Use data = sev
fit_les1 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                      sBladeArea:TidalHeight+sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                      (1|Region)+(1|MeadowId)+(1|TransectId),
                    data=sev,
                    family=Gamma(link="log"))
summary(fit_les1)
# keep Region in the random structure
### Model selection by AIC
# remove blade area : tidal height interaction
fit_les1.1 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sDensityShoots:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=Gamma(link="log"))
# remove shoot density : tidal height interaction
fit_les1.2 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sBladeArea:TidalHeight+sEpiphytePerAreaMean:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=Gamma(link="log"))
# remove epiphyte load : tidal height interaction
fit_les1.3 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        sBladeArea:TidalHeight+sDensityShoots:TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=Gamma(link="log"))
# no interactions
fit_les1.4 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+TidalHeight+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=Gamma(link="log"))
# no tidal height (or interactions)
fit_les1.5 <- glmmTMB(LesionArea~sBladeArea+sDensityShoots+sEpiphytePerAreaMean+
                        (1|Region)+(1|MeadowId)+(1|TransectId),
                      data=sev,
                      family=Gamma(link="log"))

df.AIC <- AIC(fit_les1,fit_les1.1,fit_les1.2,fit_les1.3,fit_les1.4,
              fit_les1.5)
df.AIC$deltaAIC <- df.AIC$AIC-min(df.AIC$AIC)
df.AIC$likelihood <- exp(-df.AIC$deltaAIC/2)
df.AIC$weight <- df.AIC$likelihood/sum(df.AIC$likelihood)
df.AIC
# Based on AIC and weight, use fit_les1.3, no epiphyte:tidal height interaction
summary(fit_les1.3)
performance(fit_les1.3)
# Visualize effect sizes
l_all_names <- c("Shoot density * Tidal height (upper)", "Leaf area * Tidal height (upper)",
                  "Tidal height (upper)","Epiphyte load","Shoot density","Leaf area")
l_all_plot <- plot_model(fit_les1.3,
                          axis.labels =l_all_names,
                          title="",
                          show.p = TRUE,
                          show.values = TRUE,
                          value.offset = 0.4,
                          #transform = NULL
                          #axis.lim = c(0.1,10),
                          group.terms = c(1,1,2,2,1,1))
d <- l_all_plot+theme_bw()+
  geom_hline(yintercept = 1,linetype="dashed",color="darkgrey")+
  scale_y_log10(limits=c(0.2,5))+
  scale_color_viridis_d(begin = 0, end = 0.6)+
  ylab("Scaled estimates of \nlesion area")+
  labs(tag="D")+
  theme(panel.grid = element_blank(),
        #plot.margin = margin(0,0,0,0,unit="pt"),
        axis.text = element_text(size=10))
d


# Figure of effect sizes for blade model ####
(a + b) / (c + d) / (e + f) + plot_layout(guides="collect")
ggsave(filename = "Figures/FigS6_effect_size_blade_model.jpg", width=8, height = 12)
# This is Figure S6 in the manuscript
# create high resolution version
ggsave(filename = "Figures/HighRes/FigS6_effect_size_blade_model.tiff", width=8, height = 12)
