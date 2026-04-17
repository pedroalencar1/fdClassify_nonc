
#' Multiple criteria method of FD classification
#'
#' @param vtime data frame column or vector containing \code{date} data
#' @param vtemp data frame column or vector containing daily temperature
#' @param vprec data frame column or vector containing daily precipitation
#' @param vet0 data frame column or vector containing daily ET0 (potential
#' evapotranspiration). It can be obtained with the function \code(penman_day)
#' @param veta data frame column or vector containing daily ETa (actual
#' evapotranspiration). It can be obtained with the function \code(actual_evap_day)
#' @param score a number ranging from 0 to 100. The percentile above which
#' events will be considered flash droughts.
#' @param thresholds a vector with multiple threshold values. View details for
#'  more information. We advise not to change from the default values.
#'
#' @details
#' The \code{thresholds} vector contains 8 elements respectively:
#'
#' * 1) The threshold used for cleaning data in Christian et al. method
#' (\code{Christian_clean_data_week}). It indicates the maximum value accepted
#' for the ration ETa/ET0. The default is 1.1, allowing ETa marginally higher than
#' ET0. This reduces significatly the number o missing data.
#' * 2) General threshold for anomalies. Indicates the number of standard deviations
#' to the average for the anomalies. The default value is 1.
#' * 3) SPEI threshold. The value to be considered under drought conditions. The
#' default value is the same froim Noguera et al. (2020), i.e., -1.28
#' * 4) General threshold for accumulated anomalies. Indicated the sum of weekly
#' anomalies accumulated over a period of 4 weeks. The default value is 3.
#' * 5) SESR and SPEI intensifications. The accumulated difference over 4 weeks.
#' The default value is -2.
#' * 6 and 7) EDDI intensification and recuperation. The values are the same from
#' \code{Pendergrass2020}. Default: 50 and 10 respectively
#' * 8) Delta_SESR percentile. Indicates the percentile of the difference over 4
#' weeks. We used the same threshold from Christian et al. 2020 for default (30th percentile.)
#' @md
#'
#'
#'
#'
#' @return
#' The function returns a list with two data frames. One with weekly and detailed values from the function and a second with a summary of all events identified.
#'
#' @export
#'
#' @examples
#' ET0 <- penman_day(vtime = de_tha_d$time, vwind = de_tha_d$wind_speed, vtemp = de_tha_d$temperature,
#'                   vvpd = de_tha_d$vapor_p_def, vheatflux = (de_tha_d$sensible_heat + de_tha_d$latent_heat))[,2]
#'
#' ETa <- actual_evap_day(vtime = de_tha_d$time,vlatent_heat =de_tha_d$latent_heat, vtemperature = de_tha_d$temperature)[,2]
#'
#' multi_criteria_fd <- multicriteria_fd(vtime = de_tha_d$time, vtemp = de_tha_d$temperature,
#'                                       vprec = de_tha_d$precipitation, vet0 = ET0, veta = ETa)
#'
#'

multicriteria_fd <- function(vtime, vtemp, vprec, vet0, veta, score = 0.6, d_score = 0.2,thresholds = c(1, 0, 0, 0, -2, 50, 10, 30)){


  ##### test data

  # df_d <- de_tha_d
  # score = 0.6
  # d_score = 0.2
  # thresholds = c(1, 0, 0, 0, -2, 50, 10, 30)
  #
  # ET0 <- penman_day(vtime = df_d$time, vwind = df_d$wind_speed, vtemp = df_d$temperature,
  #                   vvpd = df_d$vapor_p_def, vheatflux = (df_d$sensible_heat + df_d$latent_heat))[,2]
  # ETa <- actual_evap_day(vtime = df_d$time,vlatent_heat =df_d$latent_heat, vtemperature = df_d$temperature)[,2]
  #
  # vtime = df_d$time
  # vtemp = df_d$temperature
  # vprec = df_d$precipitation
  # vet0 = ET0
  # veta = ETa

  #build basic data frame
  df_day <- data.frame(time = vtime, temperature = vtemp, precipitation = vprec,
                       et0 = vet0, eta = veta)

  #accumulate data into weeks
  list_week_matrix <- list(NA) #list of matrixes for each variable
  df_week <- f.week(df_day[,c(1,2)],na_rm = F, f = mean, kind = 'standard')[[1]][,1]
  for (i in 2:ncol(df_day)){
    list_week_matrix[[i-1]] <- f.week(df_day[,c(1,i)],na_rm = F, f = mean, kind = 'standard')[[2]]
    df_week[,i] <- f.week(df_day[,c(1,i)],na_rm = F, f = mean, kind = 'standard')[[1]][,2]
  }
  df_week[,3:5] <- df_week[,3:5] * 7 #adjust to display sum instead of mean
  colnames(df_week) <- colnames(df_day) # rename columns
  names(list_week_matrix) <- colnames(df_day[,2:5])

  #anomaly of weekly accumulation
  df_week_anomalies <- df_week
  list_week_anomalies <- list_week_matrix
  for (i in 2:ncol(df_week)){
    list_week_anomalies[[i-1]] <- t(apply(as.data.frame(list_week_matrix[[i-1]]),1, f.anomaly))
    df_week_anomalies[,i] <- c(list_week_anomalies[[i-1]])
  }

  #join anomalies
  df_week$temp_anomaly <- df_week_anomalies$temperature
  df_week$prec_anomaly <- df_week_anomalies$precipitation
  df_week$et0_anomaly <- df_week_anomalies$et0
  df_week$eta_anomaly <- df_week_anomalies$eta

  df_week$time <- as.Date(df_week$time) #cast times from dttm to Date

  #anomaly of 4 week accumulation

  # ??

  # indexes

  ### SPEI

  df_week$deficit <- df_week$precipitation - df_week$et0 #calculate hidrological deficit
  spei_list <- f.spei(df_week$time, df_week$deficit, n = 4) #get accumulated deficit over four weeks
  df_week$spei <- spei_list[[1]]$spei
  df_week$deficit <- NULL #remove deficit column

  #EDDI
  df_week$eddi <- c(eddi_percentile(vtime = df_week$time,
                                    vet0 = df_week$et0, dist = 'tukey')
  )


  #SESR
  sesr_input_clean <- Christian_clean_data_week(vtime = df_day$time,
                                                vET0 = df_day$et0,
                                                vETa = df_day$eta,
                                                threshold = thresholds[1])

  df_week$sesr <- sesr_input_clean[[1]]$var

  n_weeks <- 4

  # accumulations and differences over intensification period
  df_week$temp_anomaly_ac <- runner::runner(df_week$temp_anomaly, k = 4, f = sum)
  df_week$prec_anomaly_ac <- runner::runner(df_week$prec_anomaly, k = 4, f = sum)
  df_week$et0_anomaly_ac <- runner::runner(df_week$et0_anomaly, k = 4, f = sum)
  df_week$eta_anomaly_ac <- runner::runner(df_week$eta_anomaly, k = 4, f = sum)
  df_week$spei_dif <- c(rep(NA, n_weeks),diff(df_week$spei, lag = n_weeks))
  df_week$eddi_dif1 <- c(rep(NA, n_weeks-2),diff(df_week$eddi, lag = n_weeks-2))
  df_week$sesr_dif <- c(rep(NA, n_weeks),diff(df_week$sesr, lag = n_weeks))

  #delta.SESR percentiles
  df_week$sesr_dif_perc <- c(t(apply(matrix(df_week$sesr_dif, nrow = 52, byrow = F),1, f.percentile)))

  #criteria

  criteria_names <- c('time', 'temp_anomaly', 'prec_anomaly', 'et0_anomaly', 'eta_anomaly',
                      'spei', 'sesr', 'eddi',
                      'prec_anomaly_ac','et0_anomaly_ac','temp_anomaly_ac',
                      'spei_dif', 'eddi_dif1', 'sesr', 'sesr_dif_perc', 'sesr_dif')


  df_criteria <- df_week[colnames(df_week) %in% criteria_names]

  df_criteria$temp_anomaly <- (df_criteria$temp_anomaly > thresholds[2])*1
  df_criteria$prec_anomaly <- (df_criteria$prec_anomaly < -thresholds[2])*1
  df_criteria$et0_anomaly <- (df_criteria$et0_anomaly > thresholds[2])*1
  df_criteria$eta_anomaly <- (df_criteria$eta_anomaly > thresholds[2])*1
  df_criteria$spei <- (df_criteria$spei < thresholds[3])*1
  df_criteria$eddi <- (df_criteria$eddi > 100*pnorm(thresholds[2],0,1))*1 #
  df_criteria$sesr <- (df_criteria$sesr < thresholds[2])*1 #

  df_criteria$temp_anomaly_ac <- (df_criteria$temp_anomaly_ac > thresholds[4])*1
  df_criteria$prec_anomaly_ac <- (df_criteria$prec_anomaly_ac < thresholds[4])*1
  df_criteria$et0_anomaly_ac <- (df_criteria$et0_anomaly_ac > thresholds[4])*1
  # df_criteria$eta_anomaly_ac <- (df_criteria$eta_anomaly_ac > thresholds[4])*1

  df_criteria$spei_dif <- (df_criteria$spei_dif < thresholds[5])*1
  df_criteria$sesr_dif <- (df_criteria$sesr_dif < thresholds[5])*1 #
  df_criteria$eddi_dif1 <- (df_criteria$eddi_dif1 < thresholds[6])*1
  df_criteria$sesr_dif_perc <- (df_criteria$sesr_dif_perc < thresholds[8])*1 #
  n_criteria <- ncol(df_criteria)-1

  # hist(df_criteria$score)
  #     #get scores and percentiles
  df_criteria$score <- rowSums(df_criteria[,2:ncol(df_criteria)], na.rm = T)/n_criteria
  df_criteria$score_percentile_global <- f.percentile(df_criteria$score)
  df_criteria$score_percentile_period <- c(t(apply(matrix(df_criteria$score,
                                                          nrow = 52, byrow = F),1,
                                                   f.percentile)))
  # ncol(df_criteria)
  # View(df_criteria)

  # classify as event
  df_events <- df_criteria[,c(1,16,17,18)]

  # df_events$is.fd <- (df_events$score_percentile_global > score)*1
  diffs <- data.frame(time = df_events$time)

  diffs$d1 <- c(0,diff(df_events$score, lag = 1))
  diffs$d2 <- c(0,0,diff(df_events$score, lag = 2))
  diffs$d3 <- c(0,0,0,diff(df_events$score, lag = 3))

  df_events$d_max <- apply(diffs[,2:4], 1, max) #get intensification
  df_events$p_max <- apply(diffs[,2:4], 1, which.max) # get intensification duration
  df_events$is.fd <- (df_events$score > score)*(df_events$d_max > d_score)

  # removes winter events
  df_events$is.fd[lubridate::month(df_events$time) %in% c(1,2,11,12)] <- 0

  #get correct length of the event (allow one week of recuperation (like christian))
  i = 1
  while (i < (nrow(df_events)-1)){
    if ((df_events$is.fd[i] == 1 & df_events$score[i+1] > score) | (df_events$is.fd[i] == 1 & df_events$score[i+2] > score)){
      df_events$is.fd[i+1] <-1
    }
    i <- i+1
  }

  # remove too short events (single week)
  for (i in 2:(nrow(df_events) - 1)){
    if (df_events$is.fd[i] == 1){
      df_events$is.fd[i] <- max(df_events$is.fd[i-1], df_events$is.fd[i+1])
    }
  }

  # add intensification period into event length
  for (i in 1:nrow(df_events)){
    if (df_events$is.fd[i] == 1){
      df_events$is.fd[(i-df_events$p_max[i]):(i-1)] <- 1
    }
  }

  index_aux3 <- rle(df_events$is.fd)$lengths
  index_aux4 <- runner::runner(index_aux3, f = sum)

  onset <- index_aux4[seq(0,length(index_aux4),2)]
  duration <- index_aux3[seq(0,length(index_aux3),2)]

  #join events that are 3 or less weeks apart

  # index_aux3 <- index_aux3[seq(1,length(index_aux3),2)]
  # index_aux4 <- index_aux4[seq(1,length(index_aux4),2)]
  #
  # if (min(index_aux3) < 4){
  #   length_short_breaks <- index_aux3[which(index_aux3 < 4)]
  #   position_short_breaks <- index_aux4[which(index_aux3 < 4)]
  #
  #   for (i in 1:length(length_short_breaks)){
  #     i = 1
  #     beg_break <- position_short_breaks[i] - length_short_breaks[i] + 1
  #     for (j in beg_break:position_short_breaks[i]){
  #       df_events$is.fd[j] <- 1
  #     }
  #   }
  # }

  # View(df_events)
  #get summary data frame
  fd.summary <- df_events[onset,]

  if (nrow(fd.summary > 0)){

    fd.summary$duration <- duration
    fd.summary$event <- 1:nrow(fd.summary)
    fd.summary$is.fd <- NULL
  }


  df_complete_data <- dplyr::left_join(df_week, df_events, by = 'time')
  df_criteria <- cbind(df_criteria, df_events$is.fd)

  output <- list(fd_timeseries = df_events, all_data = df_complete_data,
                 all_criteria = df_criteria, summary = fd.summary)

  return(output)

}


#
# ggplot(df_criteria, aes(x = 1,y = score))+
#   geom_violin()+
#   geom_jitter(shape=16, position=position_jitter(0.15))
#
# ggplot(df_criteria, aes(x = time,y = score))+
#   geom_point()
#
# ggplot(df_events, aes(x = time))+
#   geom_point(aes(y= score))+
#   geom_bar(aes(weight = is.fd))
#
#
#
