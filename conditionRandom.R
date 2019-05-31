# Conditioning random effect on the observed number of exacerbations within the previous year

imputeLastYrExac <- function (patientData){

  for (i in 1:(nrow(patientData)))

  {
    z1row <- subset (exacevents_z, ID == as.character(patientData[i, "ID"]) & Effect == "z1")
    z1 <- rnorm(1, mean = z1row$Estimate, sd = z1row$StdErrPred) #sampling from the calculated random effect dist for this individual
    z2row <- subset (exacevents_z, ID == as.character(patientData[i, "ID"]) & Effect == "z2")
    z2 <- rnorm(1, mean = z2row$Estimate, sd = z2row$StdErrPred)

    log_alpha <-   b0 +
      b_male * patientData[i, "gender"] +
      b_age10 * patientData[i, "age10"] +
      b_nowsmk * patientData[i, "nowsmk"] +
      b_oxygen * patientData[i, "oxygen"] +
      b_fev1pp100 * patientData[i, "fev1pp100"] +
      b_sgrq10 * patientData[i, "sgrq10"] +
      b_cardiovascular * patientData[i, "statin"] +
      b_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      b_LAMA * patientData[i, "LAMA"] +
      b_LABA * patientData[i, "LABA"] +
      b_ICS * patientData[i, "ICS"] +
      b_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      b_randomized_LABA * patientData[i, "randomized_LABA"] +
      b_randomized_ICS * patientData[i, "randomized_ICS"] +
      b_randomized_statin * patientData[i, "randomized_statin"] +
      b_BMI10 * patientData[i, "BMI10"]

    alpha <- exp (as.numeric(log_alpha) + z1)
    lambda <- alpha ^ gamma

    lastYrExacDraw <-  rpois(1, lambda)
    patientData[i, "imputedLastYrExacCount"] <- lastYrExacDraw

    #severity
    c_lin <-   c0 +
      c_male * patientData[i, "gender"] +
      c_age10 * patientData[i, "age10"] +
      c_nowsmk * patientData[i, "nowsmk"] +
      c_oxygen * patientData[i, "oxygen"] +
      c_fev1pp100 * patientData[i, "fev1pp100"] +
      c_sgrq10 * patientData[i, "sgrq10"] +
      c_cardiovascular * patientData[i, "statin"] +
      c_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      c_LAMA * patientData[i, "LAMA"] +
      c_LABA * patientData[i, "LABA"] +
      c_ICS * patientData[i, "ICS"] +
      c_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      c_randomized_LABA * patientData[i, "randomized_LABA"] +
      c_randomized_ICS * patientData[i, "randomized_ICS"] +
      c_randomized_statin * patientData[i, "randomized_statin"] +
      c_BMI10 * patientData[i, "BMI10"]

    OR <- exp (as.numeric(c_lin) + z2)
    pSevere <- (OR/(1+OR))

    rateSevere <- pSevere * lambda

    lastYrSevExacDraw <-  rpois(1, rateSevere)

    #making sure past hisory of hospitalization, if present, is consistent with the newly imputed history of exacerbation
    if (!is.na(patientData[i, "hosp1yr"])){
      while (TRUE) {
        if (((lastYrSevExacDraw > 0) && (patientData[i, "hosp1yr"] == 1)) || ((lastYrSevExacDraw == 0) && (patientData[i, "hosp1yr"] == 0))) {
          #message("hosp1yr is ", patientData[i, "hosp1yr"], " predicted sevCount is ", lastYrSevExacDraw)
          break
        }
        lastYrSevExacDraw <-  rpois(1, rateSevere)
      }
    }
    patientData[i, "imputedLastYrSevExacCount"] <- lastYrSevExacDraw

  }

  return (patientData)
}


validateLastYrExac <- function (patientData, loopCount = 1){
#internal validation using previous history
  for (i in 1:loopCount){
    imputedExacEvents <- imputeLastYrExac(patientData)
  }

  return(imputedExacEvents) #for debug, should be removed
}


densityLastYrExac <- function (patientData, random_sampling_N = 1e4, lastYrExacCol = "imputedLastYrExacCount", lastYrSevExacCol = "imputedLastYrSevExacCount") {

  conditionalRandEffect <- list()
  for (i in 1:(nrow(patientData)))

  {
    log_alpha <-   b0 +
      b_male * patientData[i, "gender"] +
      b_age10 * patientData[i, "age10"] +
      b_nowsmk * patientData[i, "nowsmk"] +
      b_oxygen * patientData[i, "oxygen"] +
      b_fev1pp100 * patientData[i, "fev1pp100"] +
      b_sgrq10 * patientData[i, "sgrq10"] +
      b_cardiovascular * patientData[i, "statin"] +
      b_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      b_LAMA * patientData[i, "LAMA"] +
      b_LABA * patientData[i, "LABA"] +
      b_ICS * patientData[i, "ICS"] +
      b_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      b_randomized_LABA * patientData[i, "randomized_LABA"] +
      b_randomized_ICS * patientData[i, "randomized_ICS"] +
      b_randomized_statin * patientData[i, "randomized_statin"] +
      b_BMI10 * patientData[i, "BMI10"]


    c_lin <-   c0 +
      c_male * patientData[i, "gender"] +
      c_age10 * patientData[i, "age10"] +
      c_nowsmk * patientData[i, "nowsmk"] +
      c_oxygen * patientData[i, "oxygen"] +
      c_fev1pp100 * patientData[i, "fev1pp100"] +
      c_sgrq10 * patientData[i, "sgrq10"] +
      c_cardiovascular * patientData[i, "statin"] +
      c_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      c_LAMA * patientData[i, "LAMA"] +
      c_LABA * patientData[i, "LABA"] +
      c_ICS * patientData[i, "ICS"] +
      c_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      c_randomized_LABA * patientData[i, "randomized_LABA"] +
      c_randomized_ICS * patientData[i, "randomized_ICS"] +
      c_randomized_statin * patientData[i, "randomized_statin"] +
      c_BMI10 * patientData[i, "BMI10"]

    conditionalZ <- matrix(0, nrow = random_sampling_N, ncol = 3)
    colnames(conditionalZ) <- c("weight", "z1", "z2")


    for (j in 1:random_sampling_N){

      z <- MASS::mvrnorm(1, c(0, 0), covMat)
      z1 <- z[1]
      z2 <- z[2]

      alpha <- exp (as.numeric(log_alpha) + z1)
      lambda <- as.numeric(alpha ^ gamma)

      obsLastYrExacCount <- as.numeric(patientData[i, lastYrExacCol])
      lastYrExacProb <-  dpois(obsLastYrExacCount, lambda)


      OR <- exp (as.numeric(c_lin) + z2)
      pSevere <- (OR/(1+OR))
      rateSevere <- as.numeric(pSevere * lambda)

      obsLastYrSevExacCount <- as.numeric(patientData[i, lastYrSevExacCol])
      lastYrSevExacProb <-  dpois(obsLastYrSevExacCount, rateSevere)

      conditionalZ[j, "z1"] <- z1
      conditionalZ[j, "z2"] <- z2
      conditionalZ[j, "weight"] <- lastYrExacProb*lastYrSevExacProb

    }
    ID <- as.character(patientData[i, "ID"])

    conditionalRandEffect[[ID]] <- as.data.frame(conditionalZ)

  }

  return(conditionalRandEffect)

}


validateConditionalRand <- function (patientData, conditionalZ,  random_sampling_N = 1000){


  predicted <- matrix(0, random_sampling_N, nrow(patientData))
  predicted_avg <- matrix (0, 1,  nrow(patientData))
  predicted_exac_rate <- matrix(0, random_sampling_N, nrow(patientData))
  predicted_exac_count <- matrix(0, random_sampling_N, nrow(patientData))
  predicted_severe_exac_count <- matrix(0, random_sampling_N, nrow(patientData))


  predicted_exac_probability <- matrix(0, random_sampling_N, nrow(patientData))
  predicted_severe_exac_probability <- matrix(0, random_sampling_N, nrow(patientData))

  #message('Random sampling for each individual is ', random_sampling_N)
  for (i in 1:(nrow(patientData)))

  {
    log_alpha <-   b0 +
      b_male * patientData[i, "gender"] +
      b_age10 * patientData[i, "age10"] +
      b_nowsmk * patientData[i, "nowsmk"] +
      b_oxygen * patientData[i, "oxygen"] +
      b_fev1pp100 * patientData[i, "fev1pp100"] +
      b_sgrq10 * patientData[i, "sgrq10"] +
      b_cardiovascular * patientData[i, "statin"] +
      b_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      b_LAMA * patientData[i, "LAMA"] +
      b_LABA * patientData[i, "LABA"] +
      b_ICS * patientData[i, "ICS"] +
      b_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      b_randomized_LABA * patientData[i, "randomized_LABA"] +
      b_randomized_ICS * patientData[i, "randomized_ICS"] +
      b_randomized_statin * patientData[i, "randomized_statin"] +
      b_BMI10 * patientData[i, "BMI10"]
    ID <- as.character(patientData[i, "ID"])

    z <- sample_n(conditionalZ[[ID]], random_sampling_N, replace = TRUE, weight = weight)

    #browser()
    alpha <- exp (as.numeric(log_alpha) + z[, "z1"])
    lambda <- alpha ^ gamma
    predicted_exac_rate[, i] <- lambda
    predicted_exac_count[, i] <-  as.numeric(lapply(lambda, rpois, n=1))


    patientData [i, "predicted_exac_rate"] <- mean(predicted_exac_rate[,i])
    obsRate <- patientData [i, "observedExacCount"]/patientData [i, "follow_up"]*365.25
    patientData [i, "predicted_exac_count"] <- mean(predicted_exac_count[,i])
    patientData [i, "predicted_count_low"]  <- quantile(predicted_exac_count[,i], 0.025)
    patientData [i, "predicted_count_high"] <- quantile(predicted_exac_count[, i], 0.975)

    patientData [i, "validated_exac_count"]  <- "No"
    #patientData [i, "observedExacCount"] <- (patientData[i, "observedExacCount"]!=0)
    if ((patientData [i, "observedExacCount"] <= patientData [i, "predicted_count_high"]) & (patientData [i, "observedExacCount"] >= patientData [i, "predicted_count_low"])) {
      patientData [i, "validated_exac_count"]  <- "Yes"}

    #message(summary(patientData$validated))

    #severity
    c_lin <-   c0 +
      c_male * patientData[i, "gender"] +
      c_age10 * patientData[i, "age10"] +
      c_nowsmk * patientData[i, "nowsmk"] +
      c_oxygen * patientData[i, "oxygen"] +
      c_fev1pp100 * patientData[i, "fev1pp100"] +
      c_sgrq10 * patientData[i, "sgrq10"] +
      c_cardiovascular * patientData[i, "statin"] +
      c_randomized_azithromycin * patientData[i, "randomized_azithromycin"] +
      c_LAMA * patientData[i, "LAMA"] +
      c_LABA * patientData[i, "LABA"] +
      c_ICS * patientData[i, "ICS"] +
      c_randomized_LAMA * patientData[i, "randomized_LAMA"] +
      c_randomized_LABA * patientData[i, "randomized_LABA"] +
      c_randomized_ICS * patientData[i, "randomized_ICS"] +
      c_randomized_statin * patientData[i, "randomized_statin"] +
      c_BMI10 * patientData[i, "BMI10"]

    OR <- exp (as.numeric(c_lin) + z[, "z2"])
    predicted_severe_exac_probability[, i] <- (OR/(1+OR))
    patientData [i, "predicted_severe_exac_probability"] <- mean(predicted_severe_exac_probability[,i])
    patientData [i, "predicted_severe_exac_rate"] <- patientData [i, "predicted_exac_rate"] * patientData [i, "predicted_severe_exac_probability"]

    predicted_severe_exac_count[, i] <-  as.numeric(lapply(patientData [i, "predicted_severe_exac_rate"], rpois, n=1))
    #print(predicted_severe_exac_count)
    patientData [i, "predicted_severe_exac_count"] <- mean(predicted_severe_exac_count[,i])


  }


  patientData$validated_exac_count <- as.factor (patientData$validated_exac_count)


  #  val.prob.ci.2(patientData$exac_1yr_prob, patientData$hadExacyr2, CL.smooth=TRUE, logistic.cal=TRUE, lty.log=9,
  #                col.log="red", lwd.log=1.5, col.ideal=colors()[10], lwd.ideal=0.5)


  return(patientData)

}

