
## ---------------------------
##
## Script name: admb_data_prep.R
##
## Purpose of script: Generate ADMB data inputs from simulated population
##
## Author: Claudia
##
##
##
## ---------------------------


admb_data_prep <- function(minage, maxage, fyr, lyr, fage_adult, nckmr, abund, pdeath,
                           fishery_dead, recruits, comps_option = "com_catch", 
                           prop_sampled = 1, rind_sd = 0.01, fec_known = TRUE,
                           beta = 1, comp_bias = FALSE,
                           RI_seed, CV_seed){
  
  ### survival
  pd <- pdeath %>%
    as.data.frame() %>%
    mutate(Year = 1:119) %>% filter(Year %in% (fyr-maxage+1):lyr) %>% 
    dplyr::select(-Year) %>% as.matrix()
  
  S <- as.matrix(1-pd)
  rownames(S) <- (fyr-maxage+1):lyr # year labels
  
  ### maturity
  m1 <- -9.778
  m2 <- 2.513
  
  Af = exp(m1+m2*1:maxage)/(1+exp(m1+m2*1:maxage))
  Af[1] <- 0 # age 1 not mature
  
  #### fecundity
  Fec <- 0.3 * (1:maxage)^1.4935
  FecM <- Fec # male fec same as female in the EM
  
  ### recruitment
  # set seed or RI variability
  set.seed(RI_seed)
  recr <- recruits[(fyr-maxage+1):lyr]
  recr <- recr / max(recr)
  rind_error <- stats::rnorm(n = length(recr), mean = 0, sd = rind_sd) - rind_sd^2 / 2
  recr_ind <- recr^beta * exp(rind_error)
  #index_cv <- rind_sd/recr
  
  ### composition data
  
  if(comps_option == "com_catch"){
    ### catch composition - commercial catch
    caa_com <- fishery_dead %>% filter(MortSource == "com", DeathY %in% fyr:lyr) %>%
      dplyr::select(-MortSource) %>%
      complete(DeathY, AgeLast = 1:33, fill = list(N = 0)) %>%
      spread(AgeLast, N, fill = 0) %>%
      dplyr::select(-DeathY)
    
    obs_comps <- caa_com
    
    if(prop_sampled < 1){
      obs_comps[] <- NA
      ncaa <- floor(prop_sampled * rowSums(caa_com)) # sample size as proportion of the commercial catch each year
      Tp <- length(fyr:lyr)
      
      # set seed for composition variation
      set.seed(CV_seed)
      for(d in 1:Tp){
        prob = unlist(prop.table(caa_com[d, ]))
        if(comp_bias == TRUE){
          # amplify most abundant year classes
          new_prob <- prob
          devs <- mean(prob)-prob
          new_prob <- new_prob - devs * 5
          new_prob[new_prob < 0] <- prob[new_prob < 0]
          new_prob <- new_prob/sum(new_prob)
          prob <- new_prob
        }
        obs_comps[d, ] <- as.list(stats::rmultinom(1, size = ncaa[d], prob = prob))
      }
      # convert to double
      obs_comps <- data.frame(lapply(obs_comps, as.double))
    }
  } else {
    obs_comps <- abund %>%
      filter(Year %in% fyr:lyr) %>%
      group_by(Year, Age) %>%
      summarise(N = sum(Number)) %>% ungroup() %>%
      ungroup() %>%
      complete(Year, Age = 1:33, fill = list(N = 0)) %>%
      spread(Age, N, fill = 0) %>%
      dplyr::select(-Year) %>%
      as.matrix()
    # fix comps for age 1
    obs_comps[,1] <- postdd_recruits[fyr:lyr]
    ncaa <- floor(rowSums(obs_comps))
  }
  
  # change zero to small number
  obs_comps[obs_comps == 0] <- 1e-5
  obs_comps <- apply(obs_comps, 1, function(x) x/sum(x)) # no need to t() b/c apply does it
  
  ADMB_dat <- list(fyr = fyr, lyr = lyr, fage_adult = fage_adult, fage = minage, lage = maxage,
                   fckmr = 1, nckmr = nckmr, S = t(S), Af = Af, 
                   FFec = Fec, #MFec = FecM, 
                   recr_ind = recr_ind,
                   rsd = rind_sd, comps = obs_comps, ncaa = ncaa)
  return(ADMB_dat)
}
