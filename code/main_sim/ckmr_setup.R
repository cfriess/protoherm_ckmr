

## ---------------------------
##
## Script name: ckmr_setup.R
##
## Purpose of script: to specify parameters, functional forms, and initialize dfs to run gag repro simulation
##
## Author: Claudia
##
## Date Created: 2022-04-12
##
##
## Notes:
##
##
## Some of the repro parameters were set by doing a test run for 100 years (small pop size to speed things up)
## to make sure we get a stable population
## the mortality df includes 100 years of natural mortality only for that purpose
##
## ---------------------------


set.seed(seed)

### set population size
# everything is tested with a base initial pop size of 10000, so we are just using a multiplier to scale it up
# for pop_mult 256, generating the pedigree takes about 6 hours with cmat_cmd and 30 mins with cmate; uses around 7 GB of memory
pop_size_init <- 10000 * pop_mult
stock_prob_init <-  runif(nstocks, min = 1, max = 6) # prob of being in each stock for initilization
stock_prob_init <- stock_prob_init / sum(stock_prob_init)


# ----------------------------------
#
# Annual variability
#
# ----------------------------------

### create environmental deviates for recruitment variability
# this one is now only used inside the mating function to create variability on 
# expected number of female offspring
envsd <- 0.6 # magnitude of random fluctuations (was 0.3)
x <- rnorm(nrow(Z), 0, envsd)
rho <- 0.3 # autocorrelation
mdev <- vector(length = length(x))
mdev[1] <- x[1]
for(i in 2:length(mdev)) mdev[i] <- rho * mdev[i-1] + x[i]
summary(mdev)
#plot(mdev)

### age 1 natural mortality deviation
psy <- 1/7 # probability of high or low age1 mortality in a given year
dm <- rbinom(length(x), 1, psy) # binary vector indicating what years we have decreased mortality
md <- runif(sum(dm), -0.4, -0.2) # magnitude of natural mortality reduction
im <- rbinom(length(x), 1, psy) # binary vector indicating what years we have increased mortality
mi <- runif(sum(im), 0.2, 0.4) # magnitude of natural mortality increase

xmort <- rnorm(nrow(Z), 0, 0.05) # random deviates on natural mortality; gives low levels of variability
xmort[as.logical(dm)] <- md # insert low mortality year values
xmort[as.logical(im)] <- mi # insert high mortality year values
# plot(xmort)

# ----------------------------------
#
# Initialize Population and make data objects for simulation
#
# ----------------------------------

nyears <- 119
final_y <- nyears

indiv <- makeFounders_prot(pop = pop_size_init, # size of founder pop
                           stocks = stock_prob_init, 
                           maxAge = maxAge,
                           survCurv = survship,
                           propMale = prop_male
)

# density-independent age 1M
Mage[1] <- 0.4
Mage[length(Mage)] <- Mage[length(Mage)-1]
# saveRDS(Mage, "data/intermediate/Mage.rds")

### initial juvenile abundance
# for density-dependent Age1 mortality
# ran makeFounders_prot and year 1 through mating 10 times to see what the number of ages 1-2 would be
# for setting initial value given these parameters. juv0 is ~ the mean of those 10 runs
juv0 <- 7400 # for pop_size_init 10000
juv0 <- juv0 * pop_mult # scale it up

### other needed parameters and objects

# male fecundity at age - relativized to 1

if(fec_known == 1){ 
  ### male fecundity at age same as female
  mfec_exp <- ffec_exp
} else if(fec_known == 2){
  mfec_exp <- mfec_exp_high
} else if(fec_known == 3){
  mfec_exp <- mfec_exp_low
}
male_fec <- ((1:maxAge)/maxAge)^mfec_exp

# females are allowed to switch stock
# we keep track of the original stock and reapply it after a mating season, so the stock switches are not permanent
perm_stock <- indiv %>% dplyr::select(Me, Stock) %>%
  rename(StockPerm = Stock)

# table to keep track of numbers at age by year
n <- indiv %>% group_by(AgeLast, Sex) %>%
  summarize(Number = length(unique(Me))) %>%
  ungroup() %>%
  mutate(Year = 0) %>%
  rename(Age = AgeLast) %>%
  relocate(Year, Age, Sex, Number)

# object for recruitment
recruits <- vector("numeric", nyears)
postdd_recruits <- vector("numeric", nyears)
# for deterministic Age1 probability of death (post dens-dep):
postdd_pdeath <- vector("numeric", nyears)

# table for keeping track of when males change sex
change_yr <- data.frame(Me = as.character(), ChangeYr = as.numeric())

# dfs for abundance calculated at time of mating, and fecundity over time
abund <- n[0,]
fec_df <- as.data.frame(matrix(NA, nrow = final_y, ncol = max(ages) + 1))
colnames(fec_df) <- c("Year", paste0("A",ages))
fec_df$Year <- 1:final_y

### a few more objects for assigning mortality source
FageCom <- t(apply(as.data.frame(Fsim$Fcom), 1, function(x) x * llvul))
FageRec <- t(apply(as.data.frame(Fsim$Frec), 1, function(x) x * vul))
mort_source <- array(0, dim = c(nyears, maxAge, 3), 
                     dimnames = list("Year" = 1:119, "Age" = 1:33, "Source" = c("Commercial", "Recreational", "Natural")))
mort_source[,,"Commercial"] <- FageCom / Z
mort_source[,,"Recreational"] <- FageRec / Z
mort_source[,,"Natural"] <- 1 - apply(mort_source, 1:2, sum)

archive <- make_archive() # for archiving dead
