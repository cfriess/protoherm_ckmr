

## ---------------------------
##
## Script name: sim_inputs.R
##
## Purpose of script: to specify parameters and functions used in simulations and analysis
##
## Author: Claudia
##
## Date Created: 2025-03-04
##
##
## Notes:
##
##
##
## ---------------------------

library(here)

# read in commercial and recreational fishing mortality
Fsim <- read.csv(here("data", "inputs/F.csv"))


# ----------------------------------
#
# specify simulation parameters and age schedules
#
# ----------------------------------


###############################
# life history

ages <- 1:33
maxAge <- 33

### von Bertalanffy growth parameters
linf <- 1272
k <- 0.141
t0 <- -0.331

### maturity params
m1 <- -9.778
m2 <- 2.513

### sex change params
# initial sex comp
h1 <- -7.2732033
h2 <- 0.6346547
# probability of transitioning
t1 <- 13.83
t2 <- 4.508

fec_scalar <- 0.30 # trial and error to get stable population with constant M and no environmental effects on Age0 and Age1 M
lfpar <- 4 # that scales female reproductive output by sex ratio

### function to calculate female reproductive output as a function of proportion male in the population
# perc_male here is relative to unfished perc_male
lam_fun <- function(perc_male, fem_fec, shape){ fem_fec * (1 - exp(-shape * perc_male))}

# base proportion male, given Mortality and prob of transition
# this is used to calculate the proportion male relative to unfished (lam_fun)
# obtained this number from equilibrium test runs
# it is calculated as number of age 4+ males divided by total number of age 4+ individuals in the population
base_prop_male <- 0.3 # this needs to be changed if any changes to M at age or transition probs are made

### age schedules
# Age- and stock-dependent survival (from SEDAR72, except age 1 and age 33)
# length at age
len <- linf * (1-exp(-k*(1:33-t0)))
# initial proportion male at age; used to initiate the population
prop_male = exp(h1+h2*1:maxAge)/(1+exp(h1+h2*1:maxAge))
prop_male[1:3] <- 0 # ages 1 to 3 no males
# probability of transitioning at age
prob.transition <- pnorm(q=1:maxAge, mean=t1, sd=t2)
prob.transition[1:3] <- 0 # zero percent chance of transitioning for ages 1 to 3
# maturity at age
mat_fun_male = rep(1, maxAge) # all males are mature
mat_fun_female = exp(m1+m2*1:maxAge)/(1+exp(m1+m2*1:maxAge))
mat_fun_female[1] <- 0 # age 1 not mature

# female fecundity at age
ffec_exp <- 1.4935
fec_base <- fec_scalar * (1:maxAge)^ffec_exp # resulting fec at age from this function match data

# male fecunidity at age power function exponents
mfec_exp_high <- 3.5
mfec_exp_low <- 0.1

# skipped spawning parameters for logistic fct that maps probability of of spawning to age
pspawn_steep <- 1.2
pspawn_infl <- 4.5

### natural mortality at carrying capacity (from SEDAR72, except age 0 and age 1)
# age 1 is set high because we will reduce them with reduced abundance (i.e., they are density-dependent)
# age 33 is a quasi plus group. we don't kill off everyone older than 33 but we don't want to have
# a huge pool of individuals here either so that mortality is set high
Mage <- c(2.15,0.359,0.302,0.265,0.238,0.218,0.203,0.191,0.182,0.174,
          0.167,0.162,0.158,0.154,0.15,0.147,0.145,0.143,0.141,0.139,0.138,
          0.136,0.135,0.134,0.133,0.133,0.132,0.131,0.131,0.13,0.13,0.129,0.75)

### survival and survivorship (for initiating age structure)
survival <- exp(-Mage)
survship <- vector(length = length(survival)) # starts with age 1
survship[1] <- 1
for(i in 2:length(survship)) survship[i] <- survival[i-1] * survship[i-1]
survship[length(survship)] <- survship[length(survship)]/(1-survival[length(survship)])
survship <- survship / sum(survship)


###############################
# fishery

### fishery selectivity (gamma fun)
lmax <- 600 # length at full vulnerability
b <- 200 # gamma slope
p <- 0.5 * (sqrt(lmax^2+4*b^2-lmax))

vul = (len/lmax)^(lmax/p) * exp((lmax-len)/p)
# size limit is 24", ~ 610mm. but we don't deal with 
# vulnerability to capture vs vulnerability to retention or discard mortality
# we're going to manually lower the vulnerability for ages 1 to 3 to account for 
# the size limit protection effect
vul[1:3] <- vul[1:3] * 0.25


vul <- vul / max(vul)
# need to increase older age vuln to get sex ratio sufficiently low
#vul[7:34] <- vul[6] # make flat topped for now
#plot(1:33, vul)
#saveRDS(vul, "data/intermediate/vulr.rds")

# comm LL vul
a1 <- 1.2
a0 <- 5.8
llvul <- 1 / (1+exp(-a1*(ages-a0)))
# lower age 1-4 vulnerability
# this doesn't make a huge difference since comm vuln for these 
# ages is pretty low to begin with
llvul[1:3] <- llvul[1:3] * 0.4
#saveRDS(llvul, "data/intermediate/vulc.rds")

# limits on pdeath for age 1s, just to make sure nothing wonky happens
min_a1_nat_surv <- 0.3566749 #0.5108256
max_a1_nat_surv <- 3

### F and total mortality
# make fishing mortality at age matrix
fage <- rbind(
  t(apply(as.data.frame(Fsim$Frec), 1, function(x) x * vul)) +
    t(apply(as.data.frame(Fsim$Fcom), 1, function(x) x * llvul))
)
colnames(fage) <- 1:33

### total mortality, including 100 initial years
Z <- t(apply(fage,1,function(x) x + Mage))
# convert to probability of dying (1- survival)
pdeath <- 1 - exp(-Z)


