if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
  
  saveRDS(snakemake, file = snakemake@log$snake_obj)
  
  source(".Rprofile")
}

## ---------------------------
##
## Script name: ckmr_sim.R
##
## Purpose of script: run operating model for gag CKMR project
##
## Author: Claudia & Eric
##
## Date Created: 2022-04-12
##
##
## Notes:
##
## The simulation is set up to run using snakemake. If not using snakemake,
## can specify desired OM parameters in lines 66-71 and run script by hand
##
## ---------------------------

library(tidyverse)
library(fishSim)
library(Rcpp)


options(dplyr.summarise.inform = FALSE) # silence dplyr summarize warnings

# parameterize via Snakemake if available.  Note, we don't use scen_name
# with Snakemake, as we use the output file name.
if(exists("snakemake")) {
  # get values from the snakemake object
  pop_mult <- as.numeric(snakemake@params$pop_mult)           # simulation scale (make smaller for running quick tests) - 256 is full sim level
  ssf <- as.numeric(snakemake@params$ssf)                     # female spawning site fidelity
  mating <- snakemake@params$mating                           # "random" would be the faster version of polyandry mating_sys
  eofsr <- as.integer(snakemake@params$eofsr)                 # effect of sex ratio (true = reduction in expected offspring, false = skipped spawning)
  fec_known <- snakemake@params$fec_known                     # repurposing this here. TRUE means fec proportional to body mass, false means adult fec flat
  ped_rdata <- snakemake@output$rdata
  
  # seed will be set via a hash on the output name which will be different
  # for every rep
  seed_base <- snakemake@output$rdata

  seed <- rlang::hash(seed_base) %>%
    str_replace_all("[^0-9]", "") %>%
    str_sub(end = 8) %>%
    as.integer()
  
  # print out the values
  cat("Seed is: ", seed, "\n")
  cat("pop_mult is: ", pop_mult, "\n")
  cat("ssf is: ", ssf, "\n")
  cat("mating is: ", mating, "\n")
  cat("eofsr is: ", eofsr, "\n")
  cat("fec_known is: ", fec_known, "\n")
  
} else {
  # define OM parameters
  scen_name <- "scenario1" # give it a name, for saving and loading image and data
  pop_mult <- 10 # simulation scale (make smaller for running quick tests) - 256 is full sim level, takes a while
  ssf <- 0.01 # female spawning site fidelity, number <= 1
  mating <- "random" # "cmd" or "random"; "cmd" option calls a different mating function that allows degrees of nonrandom mating
  eofsr <- 1 # sperm limitation effect (0 = no effect, 1 = reduction in expected offspring, 2 = skipped spawning)
  fec_known <- 1
  seed <- 1234 # set seed for reproducibility
  ped_rdata = paste0("data/sim_peds/ped_", scen_name, ".rds")
  dir.create(dirname(ped_rdata), showWarnings = FALSE, recursive = TRUE)
}

# further OM options
nstocks <- 40 # based on estimates of available spawning habitat
# mating system. for simplicity, using either fully random or fully (within season) monogamous
mating_sys <- ifelse(mating == "random", "polyandry", "monogamy")
rdev_by_stock <- TRUE # whether or not recruitment deviation is allowed to vary by stock
sex_ratio_by_stock <- TRUE # whether or not sex ratio is calculated by stock

### source functions
source("code/helpers/mod_funs.R")
  
### source inputs and setup scripts that specify
source("code/main_sim/sim_inputs.R")
source("code/main_sim/ckmr_setup.R")

# first and last year for population loop
start_y <- 1
end_y <- final_y
  
### for testing, can turn off fishing mortality and environmental variability like so:
# mdev[] <- 0
# xmort[] <- 0
# fage[] <- 0
# Z <- t(apply(fage,1,function(x) x + Mage))
# pdeath <- 1 - exp(-Z)

t_0 <- Sys.time()
### population loop
for(y in start_y:end_y){ 
  
  ### sex changing
  ss_res <- rcpp_sexSwitch(indiv, prob.transition, maxAge, y)
  if(y >= 70) change_yr <- rbind(change_yr, ss_res)
  
  ### mating
  mate_res <- cmate(indiv, year = y, firstBreed = 2,
                     maxClutch = Inf, maleCurve = mat_fun_male, femaleCurve = mat_fun_female,
                     nstocks = nstocks, mating = mating, eofsr = eofsr,
                     rdev_by_stock = rdev_by_stock, sex_ratio_by_stock = sex_ratio_by_stock,
                     mating_sys = mating_sys)
  offspring <- mate_res[[1]]
  indiv <- mate_res[[2]] # to track lifetime reproductive output, currently not using it
  rm(mate_res)
  # get abundance before we apply mortality and advance age to coincide with when progeny are produced
  abund <- abund %>% bind_rows(indiv %>% group_by(AgeLast, Sex) %>%
                                 summarize(Number = length(unique(Me))) %>%
                                 ungroup() %>%
                                 mutate(Year = y) %>%
                                 rename(Age = AgeLast) %>%
                                 relocate(Year, Age, Sex, Number))

  ### add permanent stock, reset Female stock
  indiv <- indiv %>% left_join(perm_stock, by = "Me")
  # for all females age > 1, set stock back to permanent stock
  indiv <- indiv %>% mutate(Stock = ifelse(Sex == "F" & AgeLast > 1, StockPerm, Stock)) %>% 
      dplyr::select(-StockPerm)
  
  ### Age 1 mortality
  # set density-dependent probability of mortality for age1 based on linear relationship for M, then convert to survival and 1-survival
  cjd <- sum(!is.na(indiv$AgeLast) & indiv$AgeLast < 3)/juv0 # current juv dens, ages 1 and 2
  # dens-dep nat survival + random variation:
  a1Mdd <- max(0 + 1.75 * cjd  + xmort[y], 0)
  #a1_nat_surv <- min(max_a1_nat_surv, max(min_a1_nat_surv, a1_nat_surv)) # restrain max and min
  # add different Age1 mortality sources to get age1 Z
  a1Z <- Mage[1] + a1Mdd + fage[(y),1] 
  a1_pdeath <- 1 - exp(-a1Z) # total prob of age1 death
  pdeath[y,1] <- a1_pdeath 
  # update mortality source likelihood for age 1
  mort_source[y,1,"Commercial"] <- FageCom[y,1] / a1Z
  mort_source[y,1,"Recreational"] <- FageRec[y,1] / a1Z
  mort_source[y,1,"Natural"] <- 1 - sum(mort_source[y,1,1:2])
  
  ### get post-dd-mortality recruits & prob death
  a1 <- abund %>% filter(Sex == "F", Age == 1, Year == y) %>% pull(Number)
  postdd_recruits[y] <- a1 * exp(-a1Mdd)
  postdd_pdeath[y] <- 1 - exp(-Mage[1]- fage[(y),1])
  # to confirm the following should be true: pdeath = 1 - densdep_surv * densindep_surv:
  ## pdeath[y,1] == 1 - exp(-a1Mdd) * exp(-Mage[1]- fage[(y),1])
  
  ### apply mortality
  rcpp_mort(indiv, y, pdeath[y,], maxAge)
  ### assign source of mortality
  newly_dead <- as.numeric(!is.na(indiv$DeathY) & indiv$DeathY == y) # 0/1 index of length nrow(indiv) where 1s died this year
  mort_probs <- mort_source[y,,] # age x mortality source matrix of relative mortality contributors for year y
  rcpp_mort_source(indiv, mort_probs, newly_dead, maxAge) # assign mortality source. 1 = com, 2 = rec, 3 = natural
  
  ### archive dead individuals 
  # take only snapshots at first to save space
  # then start tracking dead consistently in y 70 since we don't need pedigrees further back than that
  if(y %in% c(seq(10,70,10), 71:119)){
    archive <- archive_dead2(indiv = indiv, archive = archive)
  }
  indiv <- remove_dead(indiv = indiv)
  # recruits after dens-dep and dens-indep mortality
  recruits[y] <- nrow(indiv[indiv$AgeLast == 1,]) 
  
  ### increment age
  indiv <- birthdays2(indiv)
  ### add offspring - these are the new age 1s
  indiv <- rbind(indiv, offspring)
  
  ### stock switching
  # first, remake perm_stock df so we can reset stock after next year's mating
  perm_stock <- indiv %>% dplyr::select(Me, Stock) %>%
    rename(StockPerm = Stock)
  
  # shuffle females for next breeding year
  indiv <- fem_stock_switch(indiv, switch_prob = ssf)
  
  ### add to annual n at age table
  n <- n %>% bind_rows(indiv %>% group_by(AgeLast, Sex) %>%
                         summarize(Number = length(unique(Me))) %>%
                         ungroup() %>%
                         mutate(Year = y) %>%
                         rename(Age = AgeLast) %>%
                         relocate(Year, Age, Sex, Number))
  
  #cat( sprintf( '\rIteration %i complete', y))
  print(paste("it", y, ", adults = ", nrow(indiv[indiv$AgeLast > 3,]), 
              ", a1_pdeath = ", round(pdeath[y,1],2), 
              ", nage0 = ", n$Number[n$Year == y & n$Age == 1],
              #", time = ", Sys.time()
              ", jdens = ", round(cjd, 3)))
}
t_1 <- Sys.time()

archive <- rbind(archive, indiv)
indiv <- archive
rm(archive)

### data frame of all dead individuals
# that could potentially be sampled, with as many as 10 yrs of sampling
# not tracking age 1s that died from natural mortality to keep the df size manageable
recent_indiv <- indiv %>% 
  filter(DeathY >= 110)
fishery_dead <- indiv %>% 
  filter(DeathY >= 80, MortSource %in% 1:2) %>%
  group_by(AgeLast, DeathY, MortSource) %>% summarize(N = n()) %>%
  ungroup() %>%
  mutate(MortSource = ifelse(MortSource == 1, "com", "rec"))

### save results
saveRDS(
  list(
    pdeath = pdeath,
    abund_yrend = n,
    fec_df = fec_df,
    abund = abund,
    recent_indiv = recent_indiv,
    fishery_dead = fishery_dead,
    recruits = recruits,
    postdd_recruits = postdd_recruits,
    postdd_pdeath = postdd_pdeath
  ),
  file = ped_rdata
)

