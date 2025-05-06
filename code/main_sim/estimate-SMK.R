
# In order to re-run things that were run under Snakemake on the cluster
# on a local machine, you need to get the snakemake object (stored
# as an rds file in rule estimate's log file with key snake_objs) and
# store it in the variable  LOCAL_SNAKE$snake_onj, i.e., by doing
# LOCAL_SNAKE <- list()
# LOCAL_SNAKE$snake_obj <- read_rds("snake_obj.rds")
#
# Then you also need to define LOCAL_SNAKE$in_list
# to be the rds_file from rule sample_pop that is used as the
# input to rule estimate.  
#
# Setting LOCAL_SNAKE$snake_obj to the snakemake object ensures that all
# the parameter values and the seed are the same (because the
# seed is a hash off of the input file of that snakemake object.)

# LOCAL_SNAKE workflow looks something like this:
# # get the list of outputs
# LSL <- read_rds("~/Downloads/4-23-moab-post-pizza-runs.rds")
# 
# # get the one row of the estimate tibble you want (here I am just taking row 10)
# ONEROW <- LSL$estimate_snake_obj %>% slice(2801)
# 
# # source the necessary functions:
# source("code/helpers/harvest_outputs.R")
# 
# # make the LOCAL_SNAKE LIST
# LOCAL_SNAKE <- make_local_snake(ONEROW, LSL$sample_pop_out_list)
# LOCAL_SNAKE <- snake_obj_list[4,]

# For working on Eric's laptop:
# Sys.setenv(ADMB_HOME="/Users/eriq/Documents/others_code/admb-12.3")

if(exists("snakemake") && !exists("LOCAL_SNAKE")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")

  saveRDS(snakemake, file = snakemake@log$snake_obj)
  source(".Rprofile")
}

## ---------------------------
##
## Script name: estimate-SMK.R
##
## Purpose of script: run sampling and estimation model for gag CKMR project
##
## Author: Claudia & Eric
##
## Date Created: 2024-01-16
##
## Notes:
##
##
## This model has undergone several iterations. We originally incorporated uncertainty
## and bias in fisheries catch data but decided to use perfect fishery data to
## make results more easily interpretable. Some of the original inputs got
## repurposed over the evolution of the study because easier for the purpose of the snakefile paths
## to use existing variables that were no longer being used for their original purpose
##
## ---------------------------


library(tidyverse)
library(fishSim)
library(Rcpp)
library(cowplot)

# ----------------------------------
#
# Specify Estimation Model Parameters if not using snakemake
#
# ----------------------------------

### Define EM parameters

# these can be used to introduce uncertainty/ bias in fisheries data
rind_sd = 0.01 # sd for recruitment index. controls uncertainty in the recruitment index
comp_bias = FALSE # whether or not composition data are biased

# prop_sampled was originally the proportion of total commercial catch sampled for comps.
# but got repurposed to control what portion of the population we are sampling
prop_sampled = 1 # 1 (default): sample all dead individuals for ckmr, comps data based on all individuals

# beta was originally index hyperstability / hyperdepletion parameter, 
# then got repurposed for additional ckmr sampling option but was ultimately not used
beta = 1

### currently used parameters:

# fec_known controls male relative reproductive output at age
# 1 = male output equal to females, 2 = male output higher than fem, 3 = male output lower than fem
fec_known = 1 
mtHZ <- 1  # default value for the true mtDNA heterozygosity
mtHZ_est <- 1 # default value for the estimation model mtDNA

# et controls whether whether sex transition sd and male repro output are estimated in admb
# 0 = fix both, 1 = fix male repro output only, 2 = fix sex transition only, 3 = estimate both
et <- 0 

# sample_by_sex controls how CKMR sample is split by sex
# 1 = use all ages to calculate current sex ratio, then split total sample by sex ratio (default)
# 2 = use age 2+ sex ratio, 3 = no males are sampled, 4 = sampling is random by sex
sample_by_sex <- 1
                    
CKMR_seed <- 1 # each pedigree is subsampled multiple times. this controls the random seed for subsampling
RI_seed <- 1 # recruitment index seed. currently not used
CV_seed <- 100 # composition data seed. currently not used

nsampyrs <- 3 # how many years of CKMR sampling
ckmr_ssmult <- 1 # multiplier for CKMR sample size
nll_type <- 1 # whether (1; default) or not (0) to include POPs in the CKMR likelihood

seed <- 8765
set.seed(seed + CKMR_seed)

# ----------------------------------
#
# Setup
#
# ----------------------------------

admb_output_dir <- "data/admb_run_directory"

# make list to hold data objects that we will want to return.
retObjs <- list()

# read in some simulation inputs
vulr <- readRDS("data/intermediate/vulr.rds")
vulc <- readRDS("data/intermediate/vulc.rds")
Mage <- readRDS("data/intermediate/Mage.rds") 
Ffish <- read.csv("data/inputs/F.csv") # to get starting value for mean F


scen_name <- "scenario1"
in_list <- paste0("data/sim_peds/ped_", scen_name, ".rds")
outObjsFile <- paste0("data/results/EM-All-Results_", scen_name, ".rds")
outkeyObjsFile <- paste0("data/results/EM-Key_Results_", scen_name, ".rds")

if(exists("LOCAL_SNAKE")) {
  snakemake <- LOCAL_SNAKE$snake_obj[[1]]
}
if(exists("snakemake")) {
  in_list <- snakemake@input$in_list
  
  ckmr_ssmult <- as.numeric(snakemake@params$ckmr_ssmult)
  rind_sd = as.numeric(snakemake@params$rind_sd)
  prop_sampled = as.numeric(snakemake@params$prop_sampled)
  fec_known = as.integer(snakemake@params$fec_known)
  beta = as.numeric(snakemake@params$beta)
  comp_bias = as.logical(snakemake@params$comp_bias)
  mtHZ = as.numeric(snakemake@params$hzt)
  mtHZ_est = as.numeric(snakemake@params$hze)
  et = as.numeric(snakemake@params$et)
  sample_by_sex = as.integer(snakemake@params$sample_by_sex)
  CKMR_seed = as.integer(snakemake@params$ckmr_seed)
  RI_seed = as.integer(snakemake@params$ri_seed)
  CV_seed = as.integer(snakemake@params$cv_seed)
  nll_type = as.integer(snakemake@params$nll_type)
  nsampyrs = as.integer(snakemake@params$nsampyrs)
  admb_output_dir <- snakemake@output$admb_arena
  outObjsFile <- snakemake@output$rds
  outkeyObjsFile <- snakemake@output$key_obs
  
  seed_base <- in_list  # use same seed for each pedigree+sample replicate
  seed <- rlang::hash(seed_base) %>%
    str_replace_all("[^0-9]", "") %>%
    str_sub(end = 8) %>%
    as.integer()

  set.seed(seed + CKMR_seed)
  
  # print out the values
  cat("seed: ", seed, "\n")
  cat("rind_sd: ", rind_sd, "\n")
  cat("nsampyrs: ", nsampyrs, "\n")
  cat("prop_sampled is: ", prop_sampled, "\n")
  cat("fec_known is: ", fec_known, "\n")
  cat("beta is: ", beta, "\n")
  cat("comp_bias is: ", comp_bias, "\n")
  cat("et is: ", et, "\n")
  cat("mtHZ is: ", mtHZ, "\n")
  cat("mtHZ_est is: ", mtHZ_est, "\n")
  
}


if(beta == 1){ # base case
  stock_match = 1
  sample_nearshore = FALSE
} else if(beta == 2){
  stock_match = 1
  sample_nearshore = TRUE
} else if(beta == 3){
  stock_match = 0.2
  sample_nearshore = FALSE
} else {
  stock_match = 0.2
  sample_nearshore = TRUE
}


# ---------------------------------------------
#
# read in results from the simulation and transfer values to variables
#
# ---------------------------------------------


if(exists("LOCAL_SNAKE")) {
  results_list <- readRDS(paste0("C:/Users/claudia.friess/Downloads/",LOCAL_SNAKE$ped_path)) # LOCAL_SNAKE$in_list #
} else {
  results_list <- readRDS(in_list)
}

abund <- results_list$abund
recent_indiv <- results_list$recent_indiv
fishery_dead <- results_list$fishery_dead
postdd_pdeath <- results_list$postdd_pdeath

# this pdeath includes the dens-dep mortality component for age1:
pdeath <- results_list$pdeath
# this pdeath includes only the dens-independent component:
# this is only used to calculate survival which we are not using for estimation
# unless we want to fix survival
pdeath_ndd <- pdeath 
pdeath_ndd[,1] <- postdd_pdeath

# these recruits reflect post-mortality age1s:
postmort_recruits <- results_list$recruits
# these recruits reflect post-dd mortality age1s:
postdd_recruits <- results_list$postdd_recruits

# save the inlist objects to retObjs
retObjs$abund <- abund
retObjs$pdeath <- pdeath


# -----------------------------------------------
#
# Sampling Model
#
# -----------------------------------------------

source("code/helpers/admb_data_prep.R")
source("code/helpers/compile_pairs_quickin.R")
source("code/helpers/get_ckmr_samps.R")

# sample entire DEAD population or just (dead) catch for ckmr:
# my limited tests suggest this does not have much of an effect
# also set the data source for the composition data, while we're at it
if(prop_sampled == 1){
  pop <- recent_indiv # all dead
  comps_option = "full_pop" # comps based on full pop
} else {
  pop <- recent_indiv %>% filter(MortSource %in% 1:2) # only fishery dead
  comps_option = "com_catch" # comps based on commercial catch
}

### generate CKMR samples
samps <- get_ckmr_samps(pop = pop, 
                        abund = abund, 
                        n_sampyrs = nsampyrs, nsamp_mult = ckmr_ssmult, 
                        by_sex = sample_by_sex,
                        stock_match = stock_match, 
                        sample_nearshore = sample_nearshore)
#table(samps$Sex, samps$AgeLast)

qsc_pairs <- quickin(samps, max_gen = 2)

samp_summary <- samps %>% group_by(DeathY, Sex, AgeLast) %>% summarize(N = n())

### compile_pairs_quickin()
cp_out_list <- compile_pairs_quickin(samps, qsc_pairs, nll_type)

### prep CKMR data inputs into EM
kp_counts <- cp_out_list$regular %>%
  filter(c1 != c2)  # cross-cohort comparisons only

extended_kp_counts <- cp_out_list$extended %>%
  filter(c1 < c2)

# pause and look at the totals. Note that there are PO and MO pairs
# that we have not yet used in the likelihood.
extended_count_summary <- extended_kp_counts %>% 
  summarise(across(.cols = PO:Not, .fns = sum))
extended_count_summary

# can chose to save or not
# retObjs$cp_out_list <- cp_out_list
# retObjs$extended_count_summary <- extended_count_summary
# retObjs$extended_kp_counts <- extended_kp_counts
# retObjs$kp_counts <- kp_counts

#### A bunch of extra stuff for dealing with mtDNA heterozygosity  ####

# after all the other stuff got done there, we now
# are going to simulate whether each of the HSD pairs has
# the same or different mtDNA according to the heterozygosity of the
# mtDNA markers that we have.

# set a seed here so it is consistent for any 

# first, simulate how many are the same:
kpc3 <- kp_counts %>%
  mutate(
    # simulate the number of HSD pairs in each row that actually having
    # matching mtDNA because the HZ is not 1.0.
    num_matching_mtdna_hsp = rbinom(n = n(), size = HSD, prob = 1 - mtHZ),
    num_matching_mtdna_pop = rbinom(n = n(), size = PO, prob = 1 - mtHZ)
  )

# and record that:
retObjs$kpc3 <- kpc3

# then make kpc4, in which we have swapped the same HSDs to the HSS column and PO to MO
kpc4 <- kpc3 %>%
  mutate(
    HSD = HSD - num_matching_mtdna_hsp,
    HSS = HSS + num_matching_mtdna_hsp,
    PO = PO - num_matching_mtdna_pop,
    MO = MO + num_matching_mtdna_pop
  ) %>%
  select(-num_matching_mtdna_hsp, -num_matching_mtdna_pop)

# and record that
retObjs$kpc4 <- kpc4 

# -------------------------------------------
### generate other admb inputs data objects
fyr <- min(kp_counts$c1)
lyr <- max(kp_counts$c2)
minage <- 1
maxage <- 33
fage_adult <- 3 # for summing female adult abundance
recruits <- postdd_recruits # other option is postmort_recruits
retObjs$recruits <- recruits

index_weight <- comps_weight <- ckmr_weight <- 1

admb_data_list <- admb_data_prep(minage = minage, maxage = maxage,
                                 fyr = fyr, lyr = lyr, fage_adult = fage_adult,
                                 nckmr = nrow(kpc4),
                                 abund = abund, pdeath = pdeath_ndd, 
                                 fishery_dead = fishery_dead,
                                 recruits = recruits, comps_option = comps_option, 
                                 prop_sampled = prop_sampled, rind_sd = rind_sd,
                                 fec_known = fec_known, beta = 1, 
                                 comp_bias = comp_bias,
                                 RI_seed = RI_seed + seed, 
                                 CV_seed = CV_seed + seed)

ncaa <- admb_data_list$ncaa # number of age comps samples
admb_data_list$vulr <- vulr
admb_data_list$vulc <- vulc
admb_data_list$Mage <- Mage

# added estimating male relative repro output, so retooling the snakefile 'et' parameter
# to handle both sex transition and male repro output estimation flag
admb_et <- admb_ef <- -1 # base case, both are fixed
if(et == 1){ # estimate only sex transition
  admb_et <- 1 
} else if(et == 2){ # estimate only male fecundity
  admb_ef <- 1 
} else if(et == 3){ # both are estimated
  admb_et <- admb_ef <- 1 
}

admb_data_list <- c(admb_data_list[!names(admb_data_list) %in% "ncaa"], as.list(kpc4), 
                    mtHZ_est = mtHZ_est, et = admb_et, ef = admb_ef,
                    comps_weight = comps_weight, index_weight = index_weight, 
                    ckmr_weight = ckmr_weight)
admb_data_list$nll_type <- nll_type

# ballpark initial composition data weight
admb_data_list$comps_weight <-  length(admb_data_list$recr_ind) / length(admb_data_list$comps)

# reset the seed to seed here so that it is not affected by RI_seed or CV_seed
set.seed(seed)

# -----------------------------------------------
#
# Create ADMB Inputs
#
# -----------------------------------------------

library(R2admb)
source("code/helpers/admb.R") # source functions to read admb output

# create a special admb_tmp_dir so that
# multiple instances running in parallel do not overwrite one another
admb_tmp_dir <- tempfile()
dir.create(admb_tmp_dir, showWarnings = FALSE, recursive = TRUE)
cat("admb_tmp_dir: ", admb_tmp_dir, "\n")

# write ckmr.dat file
ckmrdatfile <- file.path(admb_tmp_dir, "ckmr.dat")
file.create(ckmrdatfile)
file_conn <- file(ckmrdatfile, "w")
writeLines(paste("gag.dat"), file_conn)
writeLines(paste("gag.ctl"), file_conn)
close(file_conn)

# recruitment starting values in the right neighborhood for control file
meanR <- mean(recruits[(fyr-maxage):lyr])
rball <- round(meanR, -4)
# simulated Fs
trueFrec <- mean(Ffish$Frec[(fyr-maxage):lyr])
trueFcom <- mean(Ffish$Fcom[(fyr-maxage):lyr])

# starting value for rel repro power fct exp
# note that the starting values are used if estimation for parameters is turned off
imfec_exp <- 1.4935 # base case exponent, males same as females

# when we are estimating male repro output and fec_known != 0, give it the correct starting value
if(fec_known == 2 & et %in% c(2,3) ){
  imfec_exp <- 3
} else if(fec_known == 3 & et %in% c(2,3) ){
  imfec_exp <- 0.1
}

# write control file
# these are in regular space but are log-transformed in the tpl file

admb_start_list <- list(
  irbar = rball,
  ifcom = 0.3,
  ifrec = 0.3,
  itinfl = 13.8,
  itsd = 4.5,
  mfec_exp = imfec_exp)


gagctlfile <- file.path(admb_tmp_dir, "gag.ctl")
file.create(gagctlfile)
file_conn <- file(gagctlfile, "w")
for(i in 1:length(admb_start_list)) {
  writeLines(paste("#", names(admb_start_list[i]), sep = " "), file_conn)
  writeLines(unlist(lapply(admb_start_list[i], paste, collapse=" ")), file_conn)
}
writeLines(paste("# eof"), file_conn)
writeLines("999", file_conn)

close(file_conn)

# save this
retObjs$admb_start_list <- admb_start_list

# and now we also need to copy the ckmr.tpl file to the admb_tmp_dir
tpl_file <- "ckmr_fantasy_comps.tpl"

file.copy(
  from = paste0("code/admb/", tpl_file),
  file.path(admb_tmp_dir, "ckmr.tpl"), 
  overwrite = TRUE
)

R2admb::setup_admb() # sets environment variables so admb will run from within R 
# for macs, the Sys.getenv("ADMB_HOME") varialbe is "/usr/local". for the pc, it's "c:/admb"

### while loop to run admb 1-3 times, adjusting data source weights
# the likelihood component weights are part of the data file, so need
# to re-write that in between fitting model

admb_run_cnt <- 1
admb_run_res <- data.frame(admb_run_cnt = as.numeric(rep(NA,3)),
                           index_weight = as.numeric(rep(NA,3)),
                           ckmr_weight = as.numeric(rep(NA,3)),
                           comps_weight = as.numeric(rep(NA,3)),
                           index_flag = as.logical(rep(NA,3)),
                           comps_flag = as.logical(rep(NA,3)),
                           ckmr_flag = as.logical(rep(NA,3)),
                           convergence = as.character(rep(NA,3)),
                           nll = as.numeric(rep(NA,3))
                           )

while(admb_run_cnt < 4){
  
  admb_run_res$admb_run_cnt[admb_run_cnt] <- admb_run_cnt
  admb_run_res$index_weight[admb_run_cnt] <- admb_data_list$index_weight
  admb_run_res$comps_weight[admb_run_cnt] <- admb_data_list$comps_weight
  admb_run_res$ckmr_weight[admb_run_cnt] <- admb_data_list$ckmr_weight
  
  # reset data source flags
  comps_weight_flag <- index_weight_flag <- ckmr_weight_flag <- FALSE
  
  # clean out previous run result files
  admb_files <- list.files(admb_tmp_dir)
  file.remove(setdiff(file.path(admb_tmp_dir, admb_files),
          file.path(admb_tmp_dir, c("ckmr.tpl", "gag.ctl", "ckmr.dat")))
  )
  
  # write gag.dat file
  gagdatfile <- file.path(admb_tmp_dir, "gag.dat")
  file.create(gagdatfile)
  file_conn <- file(gagdatfile, "w")
  
  
  for(i in 1:length(admb_data_list)) {
    writeLines(paste("#", names(admb_data_list[i]), sep = " "), file_conn)
    writeLines(unlist(lapply(admb_data_list[i], paste, collapse=" ")), file_conn)
  }
  writeLines(paste("# eof"), file_conn)
  writeLines("999", file_conn)
  
  close(file_conn)
  
  ### run admb
  CURDIR=getwd()
  setwd(admb_tmp_dir) # this is where our admb files live for now in the gag_ckmr project
  # I have to run this on my PC to get the admb system call to work: Sys.setenv(PATH = paste("C:/rtools40/mingw64/bin", Sys.getenv("PATH"), sep = ";"))
  system("admb ckmr", intern = TRUE) # build executable
  # capture any instances of 1 exit status
  admb_exit <- tryCatch(
    expr = {
      system("./ckmr", intern = TRUE) # run it
      0
    },
    warning = function(w){
      print(w)
      1
    }
  )
  setwd(CURDIR)
  
  ### we can have admb exit status 0 but no positive definite hessian
  # in which case there will be no .cor file
  # if admb_exit is 0, check for existance of .cor file,
  # otherwise give it exit status 
  # otherwise try upweighting index and try again
  if(!admb_exit){
    pdh <- file.exists(file.path(admb_tmp_dir, "ckmr.cor")) # check if .cor file exists
    if(!pdh) admb_exit <- 1
  }
  
  if(!admb_exit){
    ### read in report file, do data weighting
    admb_rep <- read.rep(file.path(admb_tmp_dir, "ckmr.rep"))
    admb_fit <- read.fit(file.path(admb_tmp_dir, "ckmr"))
    
    oc <- admb_rep$comps # obs age comps
    ec <- admb_rep$pahat # exp age comps
    
    oi <- admb_data_list$recr_ind # obs abundance index
    # expected recruitment = q*index
    ei <- admb_rep$q * admb_rep$rt_hat[1:(length(admb_rep$rt_hat)-1)] # exp abundance index
    
    #plot(oi)
    #lines(ei)
    
    # observed and expected ckmr pairs by subset (sampling year of the firstborn)
    ocm <- data.frame(s1 = kpc4$s1, hsd = admb_data_list$HSD) %>% 
      group_by(s1) %>% summarize(hsd = sum(hsd)) %>% pull(hsd) # obs ckmr HSD (males)
    ocf <- data.frame(s1 = kpc4$s1, hss = admb_data_list$HSS) %>%
      group_by(s1) %>% summarize(hss = sum(hss)) %>% pull(hss) # obs ckmr HSS (females)
    ocpo <- data.frame(s1 = kpc4$s1, po = admb_data_list$PO) %>%
      group_by(s1) %>% summarize(po = sum(po)) %>% pull(po)# obs ckmr PO
    ocmo <- data.frame(s1 = kpc4$s1, mo = admb_data_list$MO) %>%
      group_by(s1) %>% summarize(mo = sum(mo)) %>% pull(mo)# obs ckmr MO
    oct <- admb_data_list$HSD + admb_data_list$HSS + admb_data_list$PO + 
      admb_data_list$MO + admb_data_list$Not
    
    ecm <- data.frame(s1 = kpc4$s1, hsd = admb_rep$PHSD_vec * oct) %>% 
      group_by(s1) %>% summarize(hsd = sum(hsd)) %>% pull(hsd) # exp HSD (males)
    ecm <- pmax(.Machine$double.eps, ecm)
    ecf <- data.frame(s1 = kpc4$s1, hss = admb_rep$PHSS_vec * oct) %>%
      group_by(s1) %>% summarize(hss = sum(hss)) %>% pull(hss)  # exp HSS (females)
    ecf <- pmax(.Machine$double.eps, ecf)
    ecpo <- data.frame(s1 = kpc4$s1, po = admb_rep$PPHSD_vec * oct) %>%
      group_by(s1) %>% summarize(po = sum(po)) %>% pull(po) # exp PO
    ecpo <- pmax(.Machine$double.eps, ecpo)
    ecmo <- data.frame(s1 = kpc4$s1, mo = admb_rep$PPHSS_vec * oct) %>%
      group_by(s1) %>% summarize(mo = sum(mo)) %>% pull(mo) # exp MO
    ecmo <- pmax(.Machine$double.eps, ecmo)
    
    
    ### check for abundance data fit and upweigh, if necessary
    # used Francis 2011, T2.3
    sdnr <- sd((oi-ei) / rind_sd) # standard deviation of normalized residuals
    df <- length(oi)-1
    sdlim <- sqrt(qchisq(0.95, df)/df) # Francis 2011, p. 1133
    if(sdnr > sdlim) {
      index_weight <- index_weight * 2 # too much upweighting at once?
      index_weight_flag <- TRUE
    }
    
    ### compositions - Francis 2011 TA1.8
    oma <- apply(oc, 1, function(x) sum((1:33)*x)) # observed mean age
    ema <- apply(ec, 1, function(x) sum((1:33)*x)) # expected mean age
    emav <- apply(ec, 1, function(x) sum((1:33)^2*x)) - ema^2 # variance of expected mean age
    
    comps_weight <- 1/var((oma-ema)/sqrt(emav/ncaa))
    comps_weight_ratio <- admb_data_list$comps_weight/comps_weight
    # Francis says it typically takes a factor of at least 2 in the comps weight 
    # to affect parameter estiamtes
    if(comps_weight_ratio > 2 | comps_weight_ratio < 0.5) {
      comps_weight_flag <- TRUE
      # put some limits on how much the weight can increase or decrease
      comps_weight <- min(max(comps_weight, admb_data_list$comps_weight * 0.33), admb_data_list$comps_weight * 3)
    } else {
      comps_weight <- admb_data_list$comps_weight
    }
    
    ### ckmr - taking a mean of the male and female weights
    # using Francis 2017, Appendix A
    ckmr_resid <- (ocm-ecm)*ecm^-0.5
    ckmr_weight_m <- 1/var(ckmr_resid)
    
    ckmr_resid <- (ocf-ecf)*ecf^-0.5
    ckmr_weight_f <- 1/var(ckmr_resid)
    
    ckmr_resid <- (ocpo-ecpo)*ecpo^-0.5
    ckmr_weight_po <- 1/var(ckmr_resid)
    
    ckmr_resid <- (ocmo-ecmo)*ecmo^-0.5
    ckmr_weight_mo <- 1/var(ckmr_resid)
    
    if(nll_type){
      ckmr_weight <- mean(c(ckmr_weight_f,ckmr_weight_m, ckmr_weight_po, ckmr_weight_mo))
    } else {
      ckmr_weight <- mean(c(ckmr_weight_f,ckmr_weight_m))
    }
    
    ckmr_weight_ratio <- admb_data_list$ckmr_weight/ckmr_weight
    # using same rule of factor of 2 here as for the comps, although have not evaluated this
    # allow downweighting only
    if(ckmr_weight_ratio > 2) {
      ckmr_weight_flag <- TRUE
      # put some limits on how much the weight can increase or decrease
      ckmr_weight <- max(ckmr_weight, admb_data_list$ckmr_weight * 0.33)
    } else {
      ckmr_weight <- admb_data_list$ckmr_weight
    }
    
    admb_run_res$ckmr_flag[admb_run_cnt] <- ckmr_weight_flag
    admb_run_res$index_flag[admb_run_cnt] <- index_weight_flag
    admb_run_res$comps_flag[admb_run_cnt] <- comps_weight_flag
    admb_run_res$convergence[admb_run_cnt] <- "yes"
    admb_run_res$nll[admb_run_cnt] <- admb_fit$nlogl
    
    # check if any of the data sources need to be re-weighted
    # otherwise, set admb_run_cnt to 4 so we break out of while loop
    if(any(c(ckmr_weight_flag, comps_weight_flag, index_weight_flag))){
      admb_data_list$comps_weight <- comps_weight
      admb_data_list$index_weight <- index_weight
      admb_data_list$ckmr_weight <- ckmr_weight
      admb_run_cnt <- admb_run_cnt + 1
    } else {
      admb_run_cnt <- 4
    }
  } else {
    admb_data_list$index_weight <- admb_data_list$index_weight * 2
    admb_run_res$convergence[admb_run_cnt] <- "no"
    admb_run_cnt <- admb_run_cnt + 1
  }
}

# if more than one run that converged, re-run the one with the lowest nll
# this is here mostly to account for the fact that the last run may have been
# a non-convergence while a previous run did converge
if(sum(!is.na(admb_run_res$nll)) > 1){
  best_run <- which.min(admb_run_res$nll)
  # if best_run is the last one, we don't need to re-run
  if(best_run != max(admb_run_res$admb_run_cnt, na.rm = T)){
    admb_data_list$comps_weight <- admb_run_res$comps_weight[best_run]
    admb_data_list$index_weight <- admb_run_res$index_weight[best_run]
    admb_data_list$ckmr_weight <- admb_run_res$ckmr_weight[best_run]
    
    # clean out previous run result files
    admb_files <- list.files(admb_tmp_dir)
    file.remove(setdiff(file.path(admb_tmp_dir, admb_files),
                        file.path(admb_tmp_dir, c("ckmr.tpl", "gag.ctl", "ckmr.dat")))
    )
    
    # write gag.dat file
    gagdatfile <- file.path(admb_tmp_dir, "gag.dat")
    file.create(gagdatfile)
    file_conn <- file(gagdatfile, "w")
    
    
    for(i in 1:length(admb_data_list)) {
      writeLines(paste("#", names(admb_data_list[i]), sep = " "), file_conn)
      writeLines(unlist(lapply(admb_data_list[i], paste, collapse=" ")), file_conn)
    }
    writeLines(paste("# eof"), file_conn)
    writeLines("999", file_conn)
    
    close(file_conn)
    
    ### run admb
    CURDIR=getwd()
    setwd(admb_tmp_dir) # this is where our admb files live for now in the gag_ckmr project
    system("admb ckmr", intern = TRUE) # build executable
    # capture any instances of 1 exit status
    admb_exit <- tryCatch(
      expr = {
        system("./ckmr", intern = TRUE) # run it
        0
      },
      warning = function(w){
        print(w)
        1
      }
    )
    setwd(CURDIR)
  }
}


# save final admb_data_list
# with the starting parameters and the data list saved
# we can re-run anything later to investigate the behavior of specific runs 
retObjs$admb_data_list <- admb_data_list
retObjs$admb_run_res <- admb_run_res

rfyr <- fyr-maxage+1 # first yr the model estimates recruitment for

Nsim <- abund %>%
  filter(Age >= fage_adult) %>%
  rename(year = Year) %>%
  group_by(year, Sex) %>%
  summarise(sim_value = sum(Number)) %>%
  filter(year >= fyr, year <= lyr) %>%
  mutate(data = ifelse(Sex == "M", "Males","Females")) %>%
  dplyr::select(-Sex) %>%
  bind_rows(data.frame(year = rfyr:lyr,
                       sim_value = recruits[rfyr:lyr],
                       data = "Recruits"))

# read admb fit - only if it converged
if(!admb_exit){
  
  # read in final max gradient and estimated parameters
  admb_rep <- read.rep(file.path(admb_tmp_dir, "ckmr.rep"))
  mle_pars <- scan(file.path(admb_tmp_dir, "ckmr.par"), comment.char="#") # MLE parameter estimates
  max_grad <- as.numeric(scan(file.path(admb_tmp_dir, "ckmr.par"), what = "", comment.char = "")[16])
  
  admb_fit <- read.fit(file.path(admb_tmp_dir, "ckmr"))
  admb_std <- data.frame(par_name = admb_fit$names, value = admb_fit$est, std = admb_fit$std)
  male_ind <- which(admb_std$par_name == "Nmt")
  male_ind <- male_ind[-length(male_ind)]
  fem_ind <- which(admb_std$par_name == "ANft")
  fem_ind <- fem_ind[-length(fem_ind)]
  r_ind <- which(admb_std$par_name == "rt_hat")
  r_ind <- r_ind[-length(r_ind)]
  abund_res <- data.frame(year = c(rep(fyr:lyr, 2), (fyr-maxage+1):lyr), 
                          data = c(rep("Males", length(male_ind)),
                                  rep("Females", length(fem_ind)),
                                  rep("Recruits", length(r_ind))
                                  ),
                          est_value = c(admb_std$value[male_ind],
                                  admb_std$value[fem_ind],
                                  admb_std$value[r_ind]),
                          std = c(admb_std$std[male_ind],
                                  admb_std$std[fem_ind],
                                  admb_std$std[r_ind])) %>%
    mutate(cv = std/est_value,
           lci = est_value-1.96*std,
           uci = est_value+1.96*std) %>%
    left_join(Nsim, by = c("year","data")) %>%
    mutate(perror = (est_value - sim_value)/sim_value * 100)
  
  para_res <- admb_std[1:(which(admb_std$par_name == "wt")[1]-1),] %>%
    rename(log_value = value, log_std = std) %>%
    mutate(log_lci = log_value-1.96*log_std,
           log_uci = log_value+1.96*log_std,
           value = exp(log_value),
           lci = exp(log_lci),
           uci = exp(log_uci))
  
  if(et == 0){
    para_res$sim_value = c(meanR, trueFcom, trueFrec)
    cormat <- admb_fit$cor[1:3,1:3]
    covmat <- admb_fit$cov[1:3,1:3]
  } else if(et == 1){
    para_res$sim_value = c(meanR, trueFcom, trueFrec, 4.5)
    cormat <- admb_fit$cor[1:4,1:4]
    covmat <- admb_fit$cov[1:4,1:4]
  } else if(et == 2){
    para_res$sim_value = c(meanR, trueFcom, trueFrec, imfec_exp)
    cormat <- admb_fit$cor[1:4,1:4]
    covmat <- admb_fit$cov[1:4,1:4]
  } else if(et == 3){
    para_res$sim_value = c(meanR, trueFcom, trueFrec, 4.5, imfec_exp)
    cormat <- admb_fit$cor[1:5,1:5]
    covmat <- admb_fit$cov[1:5,1:5]
  }
  
  para_res <- para_res  %>%
    mutate(perror = (value - sim_value)/sim_value * 100)
  
} else {
  mle_pars <- max_grad <- NA
  admb_rep <- list(NA)
  abund_res <- para_res <- cormat <- covmat <- data.frame(NA)
}

# store some admb run objects:
retObjs$mle_pars <- mle_pars
retObjs$admb_rep <- admb_rep
retObjs$max_grad <- max_grad
retObjs$admb_exit_status <- admb_exit


# copy the raw files in the admb folders to the admb_output_dir
R.utils::copyDirectory(admb_tmp_dir, admb_output_dir)

sim_pop <- abund %>%
  filter(Age >= fage_adult) %>%
  group_by(Year, Sex) %>%
  summarize(Number = sum(Number)) %>%
  ungroup() %>%
  mutate(Sex = ifelse(Sex == "M", "Males","Females")) %>%
  spread(Sex, Number) %>%
  mutate(Recruits = recruits)

# true and observed CKMR data
ckmr_obs <- kpc4 %>% group_by(c2) %>%
  summarise(across(-c1, sum))
ckmr_true <- kpc3 %>% group_by(c2) %>%
  summarise(across(-c1, sum))

# retObjs is quite big and we don't need it all for the full simulation
# so we generate Key results object 
keyObjs <- list(abund_res = abund_res, admb_cor = cormat,
                admb_cov = covmat,
                para_res = para_res, sim_pop = sim_pop,
                ckmr_obs = ckmr_obs, ckmr_true = ckmr_true,
                ckmr_samps = samp_summary, admb_exit = admb_exit,
                max_grad = max_grad, admb_run_res = admb_run_res)


### save the output lists
write_rds(retObjs, file = outObjsFile, compress = "xz")
write_rds(keyObjs, file = outkeyObjsFile, compress = "xz")


