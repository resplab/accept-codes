library(tidyverse)
library(pROC)
library(haven)
library(CalibrationCurves)
source("./functions.R")
source("./prepareData.R")
source("./conditionRandom.R")


# External
message("calculating conditional random effect distributions - external")
modifiedRandExternal <- densityLastYrExac(eclipse2yr, random_sampling_N = 2e4, lastYrExacCol = "obsExac_yr1", lastYrSevExacCol = "obsExac_severe_yr1")

message("Exacerbation rates for ECLIPSE with conditional random effects - external")
ExacerbationRateResults <- exacRates(eclipse2yr, modifiedRandExternal, n_simulation=1e2, calculateROC = F)

# calibration for binary outcome
#calibratePlot(validated_ECLIPSE, outcome="binary")

# calibration for rate outcome
validated_ECLIPSE_conditionalZ <- validateConditionalRand (eclipse2yr, modifiedRandExternal, random_sampling_N = 1e3)
calibratePlot(validated_ECLIPSE_conditionalZ)

validated_ECLIPSE_conditionalZ_males <- subset(validated_ECLIPSE_conditionalZ, gender == 1)
validated_ECLIPSE_conditionalZ_females <- subset(validated_ECLIPSE_conditionalZ, gender == 0)
calibratePlot(validated_ECLIPSE_conditionalZ_males)
calibratePlot(validated_ECLIPSE_conditionalZ_females)

# Pretty ROC

library(ggthemes)
library(extrafont)
library(directlabels)
tuftefont <- choose_font(c("Gill Sans MT", "Gill Sans", "GillSans", "Verdana", "serif"), quiet = FALSE)

externalResults <- validated_ECLIPSE_conditionalZ
externalResults$hadExac <- externalResults$observedExacCount>0
externalResults$hadSevereExac <-externalResults$observedSevereExacCount>0

externalROC <- roc(predictor=externalResults$predicted_exac_rate, response=externalResults$hadExac,
               plot = T, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)

externalROCSev <-  roc(predictor=externalResults$predicted_severe_exac_rate, response=externalResults$hadSevereExac,
                   plot = F, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)


#internal
message("validating 3 trials with last year exac count imputation")
imputed <- validateLastYrExac (exacevents_wide, loopCount = 1)

message("calculating conditional random effect distributions")
modifiedRand <- densityLastYrExac(imputed, random_sampling_N = 1e4)

# message("validating 3 trials with conditional random effects")
validated_exacevents_conditionalZ <- validateConditionalRand (exacevents_wide, modifiedRand, random_sampling_N = 1e3)

message("Exacerbation rates for development set with conditional random effects - internal")
internalExacerbationRateResults <- exacRates(exacevents_wide, randomEffectDist = modifiedRand, n_simulation = 1e2)

# calibration internal

validated_exacevents_conditionalZ_males <- subset(validated_exacevents_conditionalZ, gender == 1)
validated_exacevents_conditionalZ_females <- subset(validated_exacevents_conditionalZ, gender == 0)
calibratePlot(validated_exacevents_conditionalZ_males, type = "Internal")
calibratePlot(validated_exacevents_conditionalZ_females, type = "Internal")

# pretty ROCs

internalResults <- validated_exacevents_conditionalZ
internalResults$hadExac <- internalResults$observedExacCount>0
internalResults$hadSevereExac <-internalResults$observedSevereExacCount>0

internalROC <- roc(predictor=internalResults$predicted_exac_rate, response=internalResults$hadExac,
                   plot = T, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)

internalROCSev <-  roc(predictor=internalResults$predicted_severe_exac_rate, response=internalResults$hadSevereExac,
                       plot = F, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)


rocExac <- list(external = externalROC, internal = internalROC)
rocSevExac <- list(external = externalROCSev, internal = internalROCSev)

ggroc(rocExac, size = 1) + theme_tufte(base_size = 14) + theme(legend.title = element_blank()) + xlab ("Specificity") + ylab("Sensitivity")
rocExternal <- ggroc(rocExac, size = 1.5) + theme_tufte(base_size = 30) +
  theme(legend.position = "none", axis.line = element_line(color = 'black')) + geom_abline(colour="grey", slope = 1, intercept = 1) +
  xlab ("Specificity") + ylab("Sensitivity")
ggsave("rocExac.png", plot = rocExternal, width = 10, height = 10, scale = 1, units = "in")

ggroc(rocSevExac , size = 1) + theme_tufte(base_size = 14) + theme(legend.title = element_blank()) + xlab ("Specificity") + ylab("Sensitivity")
rocExternal <- ggroc(rocSevExac , size = 1.5) + theme_tufte(base_size = 30) +
  theme(legend.position = "none", axis.line = element_line(color = 'black')) + geom_abline(colour="grey", slope = 1, intercept = 1) +
  xlab ("Specificity") + ylab("Sensitivity")
ggsave("rocSevExac.png", plot = rocExternal, width = 10, height = 10, scale = 1, units = "in")



