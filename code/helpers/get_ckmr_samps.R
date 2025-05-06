
## ---------------------------
##
## Script name: get_ckmr_samps.R
##
## Purpose of script: generate CKMR sample from the simulated population
##
## Author: Claudia
##
##
##
## ---------------------------

#'
#' @param pop FishSim-style dataframe of individuals to sample from
#' @param n_sampyrs Number of CKMR sampling years
#' @param nsamp_mult Total number of CKMR samples to collect
#' @param abund df of sex-specific abundance - used to calculate sample size as 10 * sqrt(N)
#' @param by_sex integer between 1 and 4. controls the proportion of males in the CKMR sample
#' 1 = use full sex ratio to partition CKMR N by sex, 2 = use only age 2+ females in
#' sex ratio calculation to partition CKMR N by sex, 3 = don't sample males at all,
#' 4 = sampling is random by sex
#' @param stock_match logical indicating whether or not there is sampling mismatch 
#' between stocks for juvenile and adult sampling. Defaults to 1 (no mismatch)
#' @param sample_nearshore logical controlling whether we are only sampling younger (ages 1:5)
#' individuals or not (in which case we sample ages 1 and 4+). defaults to FALSE
get_ckmr_samps <- function(pop, abund, n_sampyrs, nsamp_mult, by_sex = 1, 
                           stock_match = 1, sample_nearshore = FALSE){
  
  # calculate sample size - this is the number of samples collected each sampling year
  # number of samples is rounded to the nearest 10
  nsamp <- abund %>% 
    filter(Age >= 2, Year > max(Year) - n_sampyrs) %>%
    group_by(Year) %>%
    summarize(N = sum(Number)) %>%
    ungroup() %>%
    summarize(Nmean = mean(N)) %>%
    mutate(nsamp = floor(10 * sqrt(Nmean) / 2 / n_sampyrs / 10 ) * 10 * nsamp_mult)
  
  fage <- 1
  if(by_sex == 2) fage <- 2
  
  sr <- abund %>% 
    filter(Age >= fage) %>%
    filter(Year > max(Year) - n_sampyrs) %>%
    group_by(Year, Sex) %>%
    summarize(N = sum(Number)) %>%
    ungroup() %>%
    group_by(Sex) %>%
    summarize(Nmean = mean(N)) %>%
    spread(Sex, Nmean) %>%
    mutate(sr = M/F)
  
  samp_yrs <- (max(abund$Year) - n_sampyrs + 1):max(abund$Year)
  maleN <- ceiling(nsamp$nsamp * sr$sr)
  femaleN <- nsamp$nsamp - maleN
  nsamp_adult <- femaleN + maleN
  
  prop_stock_samp <- 0.2 + 0.8 * stock_match
  juv_stocks <- unique(pop$Stock[pop$DeathY %in% samp_yrs & pop$AgeLast == 1])
  adult_stocks <- unique(pop$Stock[pop$DeathY %in% samp_yrs & pop$AgeLast >=4])
  juv_stocks_samp <- sample(juv_stocks, floor(prop_stock_samp * length(unique(juv_stocks))))
  if(stock_match == 1){
    adult_stocks_samp <- sample(adult_stocks, floor(prop_stock_samp * length(unique(juv_stocks))))
  } else {
    stocks_matched <- floor(length(juv_stocks_samp) * stock_match)
    stocks_unmatched <- floor(prop_stock_samp * length(unique(adult_stocks))) - stocks_matched
    adult_stocks_samp_matched <- sample(intersect(juv_stocks_samp, adult_stocks), 
                                        stocks_matched )
    adult_stocks_samp_unmatched <- sample(adult_stocks[!adult_stocks %in% juv_stocks_samp], 
                                          stocks_unmatched )
    adult_stocks_samp <- c(adult_stocks_samp_matched, adult_stocks_samp_unmatched)
  }
  
  if(sample_nearshore){
    juv_age <- 1:2
    adult_age <- 3:5
  } else {
    juv_age <- 1
    adult_age <- 4:33
  }
  
  
  # sample juveniles
  ckmr_samps_juv <- get_samps(pop, n = femaleN, 
                              yrs = samp_yrs, age = juv_age, stock = juv_stocks_samp)
  # sample adults
  if(by_sex %in% 1:2){
    
    ckmr_samps_M <- get_samps(pop, n = maleN, 
                              sex = "M", yrs = samp_yrs, age = adult_age, stock = adult_stocks_samp)
    ckmr_samps_F <- get_samps(pop, n = femaleN, 
                              sex = "F", yrs = samp_yrs, age = adult_age, stock = adult_stocks_samp)
    ckmr_samps_adult <- rbind(ckmr_samps_M, ckmr_samps_F)
    
  } else if(by_sex == 3){
    ckmr_samps_adult <- get_samps(pop, n = nsamp_adult, yrs = samp_yrs,
                                  age = adult_age, sex = "F")
  } else {
    ckmr_samps_adult <- get_samps(pop, n = nsamp_adult, yrs = samp_yrs,
                                  age = adult_age)
  }
  ckmr_samps <- bind_rows(ckmr_samps_juv, ckmr_samps_adult)
  return(ckmr_samps)
}


### slightly modified capture2 function from the mod_funs.R script
# to allow sampling multiple years at once. based on fishSim code

get_samps <- function (pop2, n, yrs,  
                       sex = NULL, age = NULL, fishery = NULL,
                       mature = NULL, stock = NULL) 
{
  
  # if (!is.null(fishery)) {
  #   if(fishery == 1){
  #     is.fisherydead <- pop2[,11] %in% 1
  #   } else {
  #     is.fisherydead <- pop2[,11] %in% 2
  #   }
  # } else is.fisherydead <- pop2[,11] %in% 1:2
  
  if (!is.null(sex)) {
    is.sex <- pop2[, 2] == sex
  } else is.sex <- TRUE
  if (!is.null(age)) {
    is.age <- pop2[, 8] %in% age
  } else is.age <- TRUE
  if (!is.null(mature)) {
    stop("need to add maturity curve back in for this")
    #is.mature <- runif(nrow(pop2)) < maturityCurve[pop2[,8]]
  } else is.mature <- TRUE
  if (!is.null(stock)) {
    is.stock <- pop2[,7] %in% stock
  } else is.stock <- TRUE
  candidates <- pop2[, 6] %in% yrs & is.na(pop2[, 9]) & is.sex & is.age & is.mature & is.stock #& is.fisherydead
  # n.potential <- sum(candidates)
  # n <- floor(min(length(yrs) * n, n.potential) / length(yrs))
  if (n > 0) {
    # make sure we don't try to sample more than are there - change the
    max_prop_sampled <- 1 # change to smaller than 1 if wanting to make sure we sample less than available
    sample.loc <- do.call(c, lapply(yrs, 
                                    function(x) sample(which(pop2$DeathY[candidates] == x), 
                                                       size = min(n, (length(which(pop2$DeathY[candidates] == x)) * max_prop_sampled)) )))
    #sample.loc <- sample.int(n.potential, size = n)
    pop2[candidates, ][sample.loc, 9] <- pop2$DeathY[candidates][sample.loc]
  }
  pop2 <- pop2[!is.na(pop2$SampY), ]
  return(pop2)
}
