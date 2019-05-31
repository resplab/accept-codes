gamma	                    <- 0.9687
b0	                      <- -0.00964
b_male	                  <- -0.157
b_age10	                  <- -0.01885
b_nowsmk	                <- -0.2009
b_oxygen	                <- 0.08781
b_fev1pp100	              <- -0.4419
b_sgrq10                  <- 0.103
b_cardiovascular	        <- 0.09837
b_randomized_azithromycin <- -0.1687
b_LAMA	                  <- 0.1485
b_LABA	                  <- 0.1216
b_ICS	                    <- 0.2232
b_randomized_LAMA	        <- 0.172
b_randomized_LABA	        <- 0.1398
b_randomized_ICS	        <- -0.2452
b_randomized_statin	      <- -0.05617
b_BMI10                   <- -0.1272



c0	                      <- -3.973
c_male	                  <- 0.3889
c_age10	                  <- 0.1123
c_nowsmk                  <- 0.4025
c_oxygen                  <- 0.5558
c_fev1pp100	              <- -1.1552
c_sgrq10                  <- 0.205
c_cardiovascular	        <- 0.3255
c_randomized_azithromycin <- -0.1103
c_LAMA	                  <- -0.1385
c_LABA            	      <- 0.01246
c_ICS	                    <- 0.3879
c_randomized_LAMA	        <- 0.1074
c_randomized_LABA	        <- -0.2253
c_randomized_ICS	        <- -0.1211
c_randomized_statin	      <- 0.109
c_BMI10           	      <- -0.106


v1 	<- 0.5968
v2	<- 2.3847
cov	<- 0.147

covMat <- matrix(
  c(v1, cov, cov, v2),
  nrow = 2,
  ncol = 2
)


exacRates <- function (patientData, randomEffectDist, n_simulation = 100, calculateROC = FALSE) {

  subgroupRates <- matrix(0, ncol=18, nrow=n_simulation)
  colnames(subgroupRates) <- c("all",
                               "males",
                               "females",
                               "smokers",
                               "nonSmokers",
                               "maleSmokers",
                               "maleNonSmokers",
                               "femaleSmokers",
                               "femaleNonSmokers",
                               "severeAll",
                               "severeMales",
                               "severeFemales",
                               "severeSmokers",
                               "severeNonSmokers",
                               "severeMaleSmokers",
                               "severeMaleNonSmokers",
                               "severeFemaleSmokers",
                               "severeFemaleNonSmokers")

  results <- matrix(0, ncol=4, nrow=18)
  colnames(results) <- c("observed",
                         "predicted",
                         "predicted_lower",
                         "predicted_higher")

  rownames(results) <- c("all",
                         "males",
                         "females",
                         "smokers",
                         "nonSmokers",
                         "maleSmokers",
                         "maleNonSmokers",
                         "femaleSmokers",
                         "femaleNonSmokers",
                         "severeAll",
                         "severeMales",
                         "severeFemales",
                         "severeSmokers",
                         "severeNonSmokers",
                         "severeMaleSmokers",
                         "severeMaleNonSmokers",
                         "severeFemaleSmokers",
                         "severeFemaleNonSmokers")


  for (i in 1:n_simulation) {

    message(paste0(round(i/n_simulation*100,1), "% completed"))
    patientDataResults <- validateConditionalRand (patientData, randomEffectDist, random_sampling_N = 1e4)

    all <- patientDataResults
    males <- subset(patientDataResults, ((gender==1)))
    females <- subset(patientDataResults, ((gender==0)))
    smokers <- subset(patientDataResults, (nowsmk == 1))
    nonSmokers <- subset(patientDataResults, (nowsmk == 0))
    maleSmokers <- subset(patientDataResults, ((gender==1) & (nowsmk == 1)))
    maleNonSmokers <- subset(patientDataResults, ((gender==1) & (nowsmk == 0)))
    femaleSmokers <- subset(patientDataResults, ((gender==0) & (nowsmk == 1)))
    femaleNonSmokers <- subset(patientDataResults, ((gender==0) & (nowsmk == 0)))

    subgroupRates[i, "all"]              <- mean(all$predicted_exac_rate)
    subgroupRates[i, "males"]            <- mean(males$predicted_exac_rate)
    subgroupRates[i, "females"]          <- mean(females$predicted_exac_rate)
    subgroupRates[i, "smokers"]          <- mean(smokers$predicted_exac_rate)
    subgroupRates[i, "nonSmokers"]       <- mean(nonSmokers$predicted_exac_rate)
    subgroupRates[i, "maleSmokers"]      <- mean(maleSmokers$predicted_exac_rate)
    subgroupRates[i, "maleNonSmokers"]   <- mean(maleNonSmokers$predicted_exac_rate)
    subgroupRates[i, "femaleSmokers"]    <- mean(femaleSmokers$predicted_exac_rate)
    subgroupRates[i, "femaleNonSmokers"] <- mean(femaleNonSmokers$predicted_exac_rate)

    subgroupRates[i, "severeAll"]   <-            mean(             all$predicted_severe_exac_rate)
    subgroupRates[i, "severeMales"]   <-          mean(           males$predicted_severe_exac_rate)
    subgroupRates[i, "severeFemales"]   <-        mean(         females$predicted_severe_exac_rate)
    subgroupRates[i, "severeSmokers"]   <-        mean(         smokers$predicted_severe_exac_rate)
    subgroupRates[i, "severeNonSmokers"]   <-     mean(      nonSmokers$predicted_severe_exac_rate)
    subgroupRates[i, "severeMaleSmokers"]   <-    mean(     maleSmokers$predicted_severe_exac_rate)
    subgroupRates[i, "severeMaleNonSmokers"]   <- mean(  maleNonSmokers$predicted_severe_exac_rate)
    subgroupRates[i, "severeFemaleSmokers"]   <-  mean(   femaleSmokers$predicted_severe_exac_rate)
    subgroupRates[i, "severeFemaleNonSmokers"] <- mean(femaleNonSmokers$predicted_severe_exac_rate)

  }

  #all exacs
  results["all"            , "observed"] <- sum(all$observedExacCount)/            sum(all$follow_up)
  results["males"          , "observed"] <- sum(males$observedExacCount)/          sum(males$follow_up)
  results["females"        , "observed"] <- sum(females$observedExacCount)/        sum(females$follow_up)
  results["smokers"        , "observed"] <- sum(smokers$observedExacCount)/        sum(smokers$follow_up)
  results["nonSmokers"     , "observed"] <- sum(nonSmokers$observedExacCount)/     sum(nonSmokers$follow_up)
  results["maleSmokers"    , "observed"] <- sum(maleSmokers$observedExacCount)/    sum(maleSmokers$follow_up)
  results["maleNonSmokers" , "observed"] <- sum(maleNonSmokers$observedExacCount)/ sum(maleNonSmokers$follow_up)
  results["femaleSmokers"  , "observed"] <- sum(femaleSmokers$observedExacCount)/  sum(femaleSmokers$follow_up)
  results["femaleNonSmokers","observed"] <- sum(femaleNonSmokers$observedExacCount)/ sum(femaleNonSmokers$follow_up)

  results["all"            , "predicted"] <- mean(subgroupRates[, "all"            ])
  results["males"          , "predicted"] <- mean(subgroupRates[, "males"          ])
  results["females"        , "predicted"] <- mean(subgroupRates[, "females"        ])
  results["smokers"        , "predicted"] <- mean(subgroupRates[, "smokers"        ])
  results["nonSmokers"     , "predicted"] <- mean(subgroupRates[, "nonSmokers"     ])
  results["maleSmokers"    , "predicted"] <- mean(subgroupRates[, "maleSmokers"    ])
  results["maleNonSmokers" , "predicted"] <- mean(subgroupRates[, "maleNonSmokers" ])
  results["femaleSmokers"  , "predicted"] <- mean(subgroupRates[, "femaleSmokers"  ])
  results["femaleNonSmokers","predicted"] <- mean(subgroupRates[, "femaleNonSmokers"])

  results["all"            , "predicted_lower"] <-  quantile(subgroupRates[, "all"            ], 0.025)
  results["males"          , "predicted_lower"] <-  quantile(subgroupRates[, "males"          ], 0.025)
  results["females"        , "predicted_lower"] <-  quantile(subgroupRates[, "females"        ], 0.025)
  results["smokers"        , "predicted_lower"] <-  quantile(subgroupRates[, "smokers"        ], 0.025)
  results["nonSmokers"     , "predicted_lower"] <-  quantile(subgroupRates[, "nonSmokers"     ], 0.025)
  results["maleSmokers"    , "predicted_lower"] <-  quantile(subgroupRates[, "maleSmokers"    ], 0.025)
  results["maleNonSmokers" , "predicted_lower"] <-  quantile(subgroupRates[, "maleNonSmokers" ], 0.025)
  results["femaleSmokers"  , "predicted_lower"] <-  quantile(subgroupRates[, "femaleSmokers"  ], 0.025)
  results["femaleNonSmokers","predicted_lower"] <-  quantile(subgroupRates[, "femaleNonSmokers"],0.025)

  results["all"            , "predicted_higher"] <- quantile(subgroupRates[, "all"            ], 0.975)
  results["males"          , "predicted_higher"] <- quantile(subgroupRates[, "males"          ], 0.975)
  results["females"        , "predicted_higher"] <- quantile(subgroupRates[, "females"        ], 0.975)
  results["smokers"        , "predicted_higher"] <- quantile(subgroupRates[, "smokers"        ], 0.975)
  results["nonSmokers"     , "predicted_higher"] <- quantile(subgroupRates[, "nonSmokers"     ], 0.975)
  results["maleSmokers"    , "predicted_higher"] <- quantile(subgroupRates[, "maleSmokers"    ], 0.975)
  results["maleNonSmokers" , "predicted_higher"] <- quantile(subgroupRates[, "maleNonSmokers" ], 0.975)
  results["femaleSmokers"  , "predicted_higher"] <- quantile(subgroupRates[, "femaleSmokers"  ], 0.975)
  results["femaleNonSmokers","predicted_higher"] <- quantile(subgroupRates[, "femaleNonSmokers"],0.975)

  #severe exacs
  results["severeAll", "observed"] <-             sum(             all$observedSevereExacCount)/            sum(all$follow_up)
  results["severeMales", "observed"] <-           sum(           males$observedSevereExacCount)/          sum(males$follow_up)
  results["severeFemales", "observed"] <-         sum(         females$observedSevereExacCount)/        sum(females$follow_up)
  results["severeSmokers", "observed"] <-         sum(         smokers$observedSevereExacCount)/        sum(smokers$follow_up)
  results["severeNonSmokers", "observed"] <-      sum(      nonSmokers$observedSevereExacCount)/     sum(nonSmokers$follow_up)
  results["severeMaleSmokers", "observed"] <-     sum(     maleSmokers$observedSevereExacCount)/    sum(maleSmokers$follow_up)
  results["severeMaleNonSmokers", "observed"] <-  sum(  maleNonSmokers$observedSevereExacCount)/ sum(maleNonSmokers$follow_up)
  results["severeFemaleSmokers", "observed"] <-   sum(   femaleSmokers$observedSevereExacCount)/  sum(femaleSmokers$follow_up)
  results["severeFemaleNonSmokers","observed"] <- sum(femaleNonSmokers$observedSevereExacCount)/ sum(femaleNonSmokers$follow_up)

  results["severeAll", "predicted"] <-             mean(subgroupRates[, "severeAll"])
  results["severeMales", "predicted"] <-           mean(subgroupRates[, "severeMales"])
  results["severeFemales", "predicted"] <-         mean(subgroupRates[, "severeFemales"])
  results["severeSmokers", "predicted"] <-         mean(subgroupRates[, "severeSmokers"])
  results["severeNonSmokers", "predicted"] <-      mean(subgroupRates[, "severeNonSmokers"])
  results["severeMaleSmokers", "predicted"] <-     mean(subgroupRates[, "severeMaleSmokers"])
  results["severeMaleNonSmokers", "predicted"] <-  mean(subgroupRates[, "severeMaleNonSmokers"])
  results["severeFemaleSmokers", "predicted"] <-   mean(subgroupRates[, "severeFemaleSmokers"])
  results["severeFemaleNonSmokers","predicted"] <- mean(subgroupRates[, "severeFemaleNonSmokers"])

  results["severeAll", "predicted_lower"] <-              quantile(subgroupRates[, "severeAll"], 0.025)
  results["severeMales", "predicted_lower"] <-            quantile(subgroupRates[, "severeMales"], 0.025)
  results["severeFemales", "predicted_lower"] <-          quantile(subgroupRates[, "severeFemales"], 0.025)
  results["severeSmokers", "predicted_lower"] <-          quantile(subgroupRates[, "severeSmokers"], 0.025)
  results["severeNonSmokers", "predicted_lower"] <-       quantile(subgroupRates[, "severeNonSmokers"], 0.025)
  results["severeMaleSmokers", "predicted_lower"] <-      quantile(subgroupRates[, "severeMaleSmokers"], 0.025)
  results["severeMaleNonSmokers", "predicted_lower"] <-   quantile(subgroupRates[, "severeMaleNonSmokers"], 0.025)
  results["severeFemaleSmokers", "predicted_lower"] <-    quantile(subgroupRates[, "severeFemaleSmokers"], 0.025)
  results["severeFemaleNonSmokers","predicted_lower"] <-  quantile(subgroupRates[, "severeFemaleNonSmokers"],0.025)

  results["severeAll", "predicted_higher"] <-              quantile(subgroupRates[, "severeAll"], 0.975)
  results["severeMales", "predicted_higher"] <-            quantile(subgroupRates[, "severeMales"], 0.975)
  results["severeFemales", "predicted_higher"] <-          quantile(subgroupRates[, "severeFemales"], 0.975)
  results["severeSmokers", "predicted_higher"] <-          quantile(subgroupRates[, "severeSmokers"], 0.975)
  results["severeNonSmokers", "predicted_higher"] <-       quantile(subgroupRates[, "severeNonSmokers"], 0.975)
  results["severeMaleSmokers", "predicted_higher"] <-      quantile(subgroupRates[, "severeMaleSmokers"], 0.975)
  results["severeMaleNonSmokers", "predicted_higher"] <-   quantile(subgroupRates[, "severeMaleNonSmokers"], 0.975)
  results["severeFemaleSmokers", "predicted_higher"] <-    quantile(subgroupRates[, "severeFemaleSmokers"], 0.975)
  results["severeFemaleNonSmokers", "predicted_higher"] <- quantile(subgroupRates[, "severeFemaleNonSmokers"],0.975)

  ## ROC Curve - temporarily disabled
   patientDataResults$hadExac <- patientDataResults$observedExacCount>0
   patientDataResults$hadSevereExac <- patientDataResults$observedSevereExacCount>0

   if (calculateROC == TRUE) {
   message("ROC for exacerbations")
   roc(predictor=patientDataResults$predicted_exac_rate, response=patientDataResults$hadExac,
       plot = T, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)

   message("ROC for severe exacerbations")
   roc(predictor=patientDataResults$predicted_severe_exac_rate, response=patientDataResults$hadSevereExac,
       plot = T, ci=T, print.auc=TRUE,  boot.n=1000, ci.alpha=0.95, stratified=FALSE, show.thres=TRUE, grid=TRUE)

   }

  return(results)
}

calibratePlot <- function (patientData, outcome = "rate", type = "External"){

  if (outcome == "binary") {
      calibrationPlotData <- patientData %>% select(hadExac, predicted_exac_probability) %>%
        mutate(quantile = ntile(predicted_exac_probability, 10))

      calibrationPlot <- matrix (0, nrow=10,  ncol=3)
      colnames(calibrationPlot) <- c("decile", "mean_predicted_probability", "mean_observed_frequency")

      for (i in 1:10) {
        calibrationPlot[i, "decile"] <- i
        calibrationPlot[i, "mean_predicted_probability"] <- mean(subset(calibrationPlotData, quantile==i)$predicted_exac_probability)
        calibrationPlot[i, "mean_observed_frequency"] <- mean(subset(calibrationPlotData, quantile==i)$hadExac)
      }
      calibrationPlot <- as.data.frame(calibrationPlot)
      p <- ggplot(calibrationPlot) +
        geom_point(aes(y = mean_observed_frequency, x = mean_predicted_probability), shape="triangle") +
        geom_smooth(aes(y = mean_observed_frequency, x = mean_predicted_probability), method = loess, se = FALSE, linetype=3) +
        #xlim(0, 1) + ylim(0, 1) +
        geom_abline(linetype=2) + theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_point(data = calibrationPlotData, aes(x=predicted_exac_probability, y=as.integer(hadExac)), shape=1, size = 1) +
        geom_jitter(data = calibrationPlotData, aes(x=predicted_exac_probability, y=as.integer(hadExac)),
                    shape=1, size = 1,  width = 0, height=0.03)
      plot(p)

      val.prob.ci.2(patientData$predicted_exac_probability, patientData$hadExac, CL.smooth=TRUE, logistic.cal=TRUE, lty.log=9,
                    col.log="red", lwd.log=1.5, col.ideal=colors()[10], lwd.ideal=0.5)
  }

  if (outcome == "rate") {
    calibrationPlotData <- patientData %>% select(observedExacCount, predicted_exac_count) %>%
      mutate(quantile = ntile(predicted_exac_count, 10))

    calibrationPlot <- matrix (0, nrow=10,  ncol=5)
    colnames(calibrationPlot) <- c("decile", "mean_predicted_count", "mean_observed_count", "SE_predicted_count", "SE_observed_count")

    for (i in 1:10) {
      calibrationPlot[i, "decile"] <- i
      calibrationPlot[i, "mean_predicted_count"] <- mean(subset(calibrationPlotData, quantile==i)$predicted_exac_count)
      calibrationPlot[i, "mean_observed_count"] <- mean(subset(calibrationPlotData, quantile==i)$observedExacCount)
      calibrationPlot[i, "SE_predicted_count"] <- sd(subset(calibrationPlotData, quantile==i)$predicted_exac_count) / sqrt(length(subset(calibrationPlotData, quantile==i)$predicted_exac_count))
      calibrationPlot[i, "SE_observed_count"] <- sd(subset(calibrationPlotData, quantile==i)$observedExacCount) / sqrt(length(subset(calibrationPlotData, quantile==i)$observedExacCount))
    }
    calibrationPlot <- as.data.frame(calibrationPlot)
    plotTitle <- paste0(type, " Validation, n=", dim(patientData)[1] )
    p <- ggplot(calibrationPlot) +
      geom_point(aes(y = mean_observed_count, x = mean_predicted_count), shape="triangle") +
      geom_smooth(aes(y = mean_observed_count, x = mean_predicted_count), method = loess, se = T, linetype=2) +
      #xlim(0, 1) + ylim(0, 1) +
      geom_abline(linetype=3) + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ggtitle(plotTitle) + xlab("Predicted Rate of Exacerbation") + ylab("Observed Rate of Exacerbation")

    plot(p)

    #alternative formatting for calibration plot

    p2 <- ggplot(calibrationPlot, aes(y = mean_observed_count, x = mean_predicted_count)) +
      geom_point() +
      geom_errorbar(aes(ymin = (mean_observed_count-1.96*SE_observed_count), ymax = mean_observed_count+1.96*SE_observed_count), width = 0.2) +
      #geom_smooth(aes(y = mean_observed_count, x = mean_predicted_count), method = loess, se = T, linetype=2) +
      #xlim(0, 5) + ylim(0, 5) +
      geom_abline(linetype=3) + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ggtitle(plotTitle) + xlab("Predicted Rate of Exacerbation") + ylab("Observed Rate of Exacerbation") + theme_tufte() +
      theme(axis.line = element_line(color = 'black'))

    plot(p2)

    # Severe Exacerbation only:
    calibrationPlotData <- patientData %>% select(observedSevereExacCount, predicted_severe_exac_count) %>%
      mutate(quantile = ntile(predicted_severe_exac_count, 10))

    calibrationPlot <- matrix (0, nrow=10,  ncol=5)
    colnames(calibrationPlot) <- c("decile", "mean_predicted_count", "mean_observed_count", "SE_predicted_count", "SE_observed_count")

    for (i in 1:10) {
      calibrationPlot[i, "decile"] <- i
      calibrationPlot[i, "mean_predicted_count"] <- mean(subset(calibrationPlotData, quantile==i)$predicted_severe_exac_count)
      calibrationPlot[i, "mean_observed_count"] <- mean(subset(calibrationPlotData, quantile==i)$observedSevereExacCount)
      calibrationPlot[i, "SE_predicted_count"] <- sd(subset(calibrationPlotData, quantile==i)$predicted_exac_count) / sqrt(length(subset(calibrationPlotData, quantile==i)$predicted_exac_count))
      calibrationPlot[i, "SE_observed_count"] <- sd(subset(calibrationPlotData, quantile==i)$observedExacCount) / sqrt(length(subset(calibrationPlotData, quantile==i)$observedExacCount))
    }
    calibrationPlot <- as.data.frame(calibrationPlot)
    plotTitle <- paste0(type, " Validation, n=", dim(patientData)[1] )
    pSevere <- ggplot(calibrationPlot) +
      geom_point(aes(y = mean_observed_count, x = mean_predicted_count), shape="triangle") +
      geom_smooth(aes(y = mean_observed_count, x = mean_predicted_count), method = loess, se = T, linetype=2) +
      #xlim(0, 1) + ylim(0, 1) +
      geom_abline(linetype=3) + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ggtitle(plotTitle) + xlab("Predicted Rate of Severe Exacerbation") + ylab("Observed Rate of Severe Exacerbation")

    plot(pSevere)

    #alternative formatting for calibration plot

    p2Severe <- ggplot(calibrationPlot, aes(y = mean_observed_count, x = mean_predicted_count)) +
      geom_point() +
      geom_errorbar(aes(ymin = (mean_observed_count-1.96*SE_observed_count), ymax = mean_observed_count+1.96*SE_observed_count), width = 0.2) +
      #geom_smooth(aes(y = mean_observed_count, x = mean_predicted_count), method = loess, se = T, linetype=2) +
      #xlim(0, 5) + ylim(0, 5) +
      geom_abline(linetype=3) + theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      ggtitle(plotTitle) + xlab("Predicted Rate of Exacerbation") + ylab("Observed Rate of Exacerbation")

    plot(p2Severe)

  }
}

