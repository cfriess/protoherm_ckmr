
## ---------------------------
##
## Script name: mod_funs.R
##
## Purpose of script: Functions needed to run simulation. A number of these are modified versions 
## of the fishSim package functions while others are new functions or Rcpp versions of functions
##
## Author: Claudia
##
##
##
## ---------------------------



# --------------------------------------------------------------
# 
# Modified FishSim functions
#
# --------------------------------------------------------------

# ---------------------------
# makeFounders

# set up founder population for protogynous hermaphrodites
# assigns sex via binomial draws based on proportion male at age in the population

#' @param propMale vector of equilibrium proportion male at age
makeFounders_prot <- function (pop = 1000, 
                               stocks = c(0.3, 0.3, 
                                          0.4), maxAge = 20, 
                               survCurv = 0.7^(1:maxAge)/sum(0.7^(1:maxAge)),
                               propMale = herm_fun # equilibrium proportion male 
) 
{
  if (sum(stocks) != 1) 
    warning("stocks do not sum to 1")
  if (sum(survCurv) != 1) 
    warning("survCurv does not sum to 1")
  if (length(survCurv) != maxAge) 
    warning("survCurv and maxAge imply different maximum ages")
  indiv <- data.frame(Me = character(pop), Sex = character(pop), 
                      Dad = character(pop), Mum = character(pop), BirthY = integer(pop), 
                      DeathY = integer(pop), Stock = integer(pop), AgeLast = integer(pop), 
                      SampY = integer(pop), Cum_os = integer(pop), MortSource = integer(pop))
  indiv[, 1] <- ids::uuid(n = pop, drop_hyphens = TRUE)
  indiv[, 3] <- c(rep("founder", pop))
  indiv[, 4] <- c(rep("founder", pop))
  indiv[, 6] <- as.integer(c(rep(NA, pop)))
  indiv[, 7] <- as.integer(sample(1:length(stocks), pop, TRUE, 
                                  prob = stocks))
  indiv[, 8] <- sample.int(maxAge, pop, TRUE, prob = survCurv)
  indiv[, 5] <- 1L - indiv[, 8]
  indiv$Sex <- apply(as.data.frame(indiv$AgeLast), 1, function(x) rbinom(1,1,propMale[x]))
  indiv$Sex <- ifelse(indiv$Sex == 0, "F", "M")
  indiv[, 9] <- as.integer(c(rep(NA, pop)))
  indiv[, 10] <- rep(0, pop)
  indiv[, 11] <- as.integer(c(rep(NA, pop)))
  return(indiv)
}


# ---------------------------
# rcpp_sexSwitch

#' rcpp version of protogynous sex switch function. faster
#' This function returns only a dataframe of the ids and sex switch year 
#' for the individuals who switch sex.
#' The 'indiv' dataframe is updated without having to return it
#' because Rcpp uses shallow copies by default. So sex for the input df will automatically
#' be updated without having to return it

#' @param N individiuals data frame
#' @param tprob numeric vector, annual probability of transitioning to male at age
#' @param maxAge integer. age of the plus group
#' @param year integer. current year
cppFunction('DataFrame rcpp_sexSwitch(DataFrame N, NumericVector tprob, int maxAge, int year) {
  
  IntegerVector Age = N["AgeLast"];
  CharacterVector Sex = N["Sex"];
  CharacterVector Me = N["Me"];
  IntegerVector SexInt (Sex.size());
  IntegerVector ChangeYr (Sex.size(), 0);
 
  // SexInt is 0 for females and 1 for males
  for(int i = 0; i < Sex.size();i++) {
    if(Sex[i] == "F") {
      SexInt[i] = 0;
    }
    else {
      SexInt[i] = 1;
    }
  }
  // now subset only for the females
  IntegerVector AgeF = Age[SexInt == 0];
  CharacterVector SexF = Sex[SexInt == 0];
  IntegerVector ChangeYrF = ChangeYr[SexInt == 0];
  
  // now loop over ages and use binomial draws to decide who switches sex
  for (int a = 1; a <= maxAge; a++) {
    IntegerVector AgeSub = AgeF[AgeF == a];
    CharacterVector SexSub = SexF[AgeF == a];
    IntegerVector ChangeYrSub = ChangeYrF[AgeF == a];
    int n = AgeSub.size();
    NumericVector prob = rbinom(n, 1, tprob[a-1]);
    SexSub[prob == 1] = "M";
    ChangeYrSub[prob == 1] = year;
    SexF[AgeF == a] = SexSub;
    ChangeYrF[AgeF == a] = ChangeYrSub;
  }
  Sex[SexInt == 0] = SexF;
  ChangeYr[SexInt == 0] = ChangeYrF;
 
  DataFrame result = DataFrame::create(Named("Me") = Me[ChangeYr != 0],
                                          Named("ChangeYr") = ChangeYr[ChangeYr != 0]);
  
  return result;
  
}')


# ---------------------------
# cmate

#' gets inputs ready for rcpp_mating functions, then calls rcpp function depending on desired mating scheme
#' female reproductive output is scaled by relative male abundance to simulate the effects of sperm limitation
#' @param indiv data frame of individuals used in fishSim
#' @param year current simulation year
#' @param param firstBreed age at first breeding
#' @param maxClutch maximum number of clutch size per female
#' @param maleCurve male maturity at age vector (essentially useless here since all males are mature)
#' @param femaleCurve female maturity at age vector
#' @param mating which mating system do we want? One of "random" (random mating within stock) or "cmd" (father selection according to the compound multinomial dirichlet)
#' @param rdev_by_stock logical. does inter-annual recruitment success vary by stock?
#' @param sex_ratio_by_stock logical. is sex ratio computed for the population as a whole or for each stock separately?
#' @param mating_sys either "monogamy", "polyandry" or "mix". this determines lambda for the CMD
#' @param eofsr integer, Either 0, 1, or 2. Sperm limitation effect. 0 = no effect (same number of eggs per female produced, no skipped spawning), 1 =  expected reproductive output reduced proportionally across female ages, 3 = no effect on expected number of offspring produced per female but increases in probability of skipped spawning

cmate <- function (indiv = makeFounders(), year = "-1", firstBreed = 1,  
                   maxClutch = Inf, maleCurve, femaleCurve, nstocks, eofsr = 1,
                   mating = c("random","cmd"), rdev_by_stock = TRUE, sex_ratio_by_stock = FALSE,
                   mating_sys) {
  
  indiv$id <- 1:nrow(indiv)
  ### select mothers and fathers
  mothers <- subset(indiv, indiv[, 2] == "F" & indiv[, 8] >= firstBreed & is.na(indiv[, 6]))
  fathers <- subset(indiv, indiv[, 2] == "M" & indiv[, 8] >= firstBreed & is.na(indiv[, 6]))
  # randomly draw mothers based on age and maturity at age AgeLast
  mothers <- mothers[runif(nrow(mothers)) < femaleCurve[mothers[,8]], , drop = FALSE]
  # randomly draw fathers based on age and maturity at age AgeLast
  fathers <- fathers[runif(nrow(fathers)) < maleCurve[fathers[, 8]], , drop = FALSE]
  
  ### calculate sex ratio, expected number of female offspring, and mother clutch size
  # fecundity_at_age as a function of sex ratio (fewer males -> fewer offspring)
  # and random variation (mdev) to get variability in recruitment - either by stock or by pop as whole
  # eofsr determines whether sex ratio affects expected number of offspring or probability of a female skipping spawning
  
  # use sigmoid function to map pskip to age
  pspawn_sigmoid <- 1/(1+exp(-pspawn_steep*(ages-pspawn_infl)))
  
  if(sex_ratio_by_stock == FALSE){
    
    sex_ratio <- nrow(fathers) / (nrow(fathers) + nrow(mothers))
    
    if(rdev_by_stock){
      mdev_stock <- runif(nstocks)
      mdev_stock <- mdev_stock * exp(mdev[y]) / sum(mdev_stock) * nstocks
      if(eofsr == 0){
        current_lam <- data.frame(exp_os = fec_base) # nages x 1 df
        female_fec <- apply(current_lam, 1, function(x) x * mdev_stock) # nstocks x nages
      } else if(eofsr == 1) {
        current_lam <- data.frame(exp_os = lam_fun(sex_ratio/base_prop_male, fec_base, lfpar)) # nages x 1 df
        female_fec <- apply(current_lam, 1, function(x) x * mdev_stock) # nstocks x nages
      } else {
        current_lam <- data.frame(fec_base) # nages x 1 df
        female_fec <- apply(current_lam, 1, function(x) x * mdev_stock) # nstocks x nages
        # probability of skipping spawning
        if (eofsr == 2){ # age invariant
          pskip <- 1 - lam_fun(sex_ratio/base_prop_male, 1, lfpar)
          mothers <- mothers[runif(nrow(mothers)) > pskip, , drop = FALSE]
        } else { # skipping spawning prob higher for younger individuals
          pspawn_base <- lam_fun(sex_ratio/base_prop_male, 1, lfpar*2)
          # higher prob of spawning for younger females
          pspawn <- pspawn_sigmoid * pspawn_base
          mothers <- mothers[runif(nrow(mothers)) < pspawn[mothers[,8]], , drop = FALSE]
        }
      }
      clutch <- mapply(function(s, a) rpois(n = 1, lambda = female_fec[s, a]), mothers$Stock, mothers$AgeLast)
    } else {
      if(eofsr == 0){
        female_fec <- fec_base * exp(mdev[y])
      } else if(eofsr == 1){
        female_fec <- lam_fun(sex_ratio/base_prop_male, fec_base, lfpar) * exp(mdev[y]) # single vector length nages
      } else {
        female_fec <- fec_base * exp(mdev[y])
        # probability of skipping spawning
        if(eofsr == 2){
          pskip <- 1 - lam_fun(sex_ratio/base_prop_male, 1, lfpar)
          mothers <- mothers[runif(nrow(mothers)) > pskip, , drop = FALSE]
        } else {
          pspawn_base <- lam_fun(sex_ratio/base_prop_male, 1, lfpar*2)
          # higher prob of spawning for younger females
          pspawn <- pspawn_sigmoid * pspawn_base
          mothers <- mothers[runif(nrow(mothers)) < pspawn[mothers[,8]], , drop = FALSE]
        }
      }
      clutch <- apply(mothers, 1, function(x) rpois(n = 1, lambda = female_fec[as.numeric(x[8])]))
    }
  } else {
    
    sex_ratio <- mothers %>% group_by(Stock) %>% summarize(F = n()) %>% 
      left_join(fathers %>% group_by(Stock) %>% summarize(M = n()), by = "Stock") %>%
      complete(Stock = 1:40, fill = list(F = 0, M = 0)) %>%
      mutate(M = ifelse(is.na(M), 0, M)) %>%
      mutate(current = M/(F+M), 
             base = base_prop_male)
    
    if(eofsr == 0){
      current_lam <- matrix(rep(fec_base, nstocks), nrow = nstocks, byrow = T)
    } else if(eofsr == 1){
      current_lam <- t(apply(as.data.frame(sex_ratio), 1, function(x) lam_fun(x[4]/x[5], fec_base, lfpar))) # nstocks x nages df
    } else {
      current_lam <- matrix(rep(fec_base, nstocks), nrow = nstocks, byrow = T)
      # different probability of skipping spawning for each stock
      if(eofsr == 2){
        pskip <- 1 - lam_fun(sex_ratio$current/sex_ratio$base, 1, lfpar)
        mothers$id <- 1:nrow(mothers)
        mothers <- lapply(1:nstocks, function(S){
          tmp = mothers[mothers$Stock == S,]
          res = tmp[runif(nrow(tmp)) > pskip[S], , drop = FALSE]
        }) %>% bind_rows()
        mothers <- mothers[order(mothers$id),]
        mothers <- mothers[, -ncol(mothers)]
      } else {
        pspawn_base <- lam_fun(sex_ratio$current/sex_ratio$base, 1, lfpar*2)
        # map to nstock x nages matrix
        pspawn <- t(apply(as.data.frame(pspawn_base), 1, function(x) x * pspawn_sigmoid))
        mothers <- lapply(1:nstocks, function(S){
          tmp = mothers[mothers$Stock == S,]
          res = tmp[runif(nrow(tmp)) < pspawn[S, tmp[,8]], , drop = FALSE]
        }) %>% bind_rows()
        mothers <- mothers[order(mothers$id),]
        mothers <- mothers[, -ncol(mothers)]
      }
    }
    if(rdev_by_stock){
      mdev_stock <- runif(nstocks)
      mdev_stock <- mdev_stock * exp(mdev[y]) / sum(mdev_stock) * nstocks
      female_fec <- apply(current_lam, 2, function(x) x * mdev_stock) # nsstock x nages
    } else {
      female_fec <- current_lam * exp(mdev[y])
    }
    clutch <- mapply(function(s, a) rpois(n = 1, lambda = female_fec[s, a]), mothers$Stock, mothers$AgeLast)
  }
  
  mothers <- subset(mothers, clutch > 0) # select only mothers with more than 0 offspring at age 0
  clutch <- clutch[clutch > 0]
  clutch[clutch > maxClutch] <- maxClutch
  mothers$Cum_os <- clutch
  # probability vector of length nfathers for each male to father offspring
  probs <- male_fec[fathers$AgeLast]
  
  if(mating == "random"){
    offspring <- rcpp_mating(mothers, fathers, unname(clutch), probs, nstocks, y)
    offspring$Me <- ids::uuid(n = nrow(offspring), drop_hyphens = TRUE)
    offspring$Cum_os <- as.integer(0)
    offspring$MortSource <- as.integer(NA)
  } else if (mating == "cmd") {
    result <- rcpp_mate_cmd(nrow(mothers), nrow(fathers), clutch, probs, 
                            mothers$Stock, fathers$Stock, mothers$AgeLast, max(mothers$Stock), mating_sys)
    offspring <- data.frame(Me = ids::uuid(n = nrow(result), drop_hyphens = TRUE),
                            Sex = "F",
                            Dad = fathers$Me[result$Dad],
                            Mum = mothers$Me[result$Mum],
                            BirthY = as.integer(y),
                            DeathY = as.integer(NA),
                            Stock = as.integer(result$Stock),
                            AgeLast = as.integer(1),
                            SampY = as.integer(NA), 
                            Cum_os = as.integer(0),
                            MortSource = as.integer(NA))
    
  } else {
    stop("no such mating function")
  }
  fathers <- fathers %>% dplyr::select(-Cum_os) %>% left_join(offspring %>% 
                                                                group_by(Dad) %>% summarize(Cum_os = n()) %>% 
                                                                rename(Me = Dad), by = "Me")
  fathers$Cum_os[is.na(fathers$Cum_os)] <- 0
  indiv$Cum_os[indiv$id %in% mothers$id] <- indiv$Cum_os[indiv$id %in% mothers$id] + mothers$Cum_os
  indiv$Cum_os[indiv$id %in% fathers$id] <- indiv$Cum_os[indiv$id %in% fathers$id] + fathers$Cum_os
  
  return(list(offspring = offspring, indiv = indiv[,-ncol(indiv)]))
}


# ---------------------------
# rcpp_mating

# simple Rcpp version of random mating function, called from the cmate R function

#' @param mothers data frame females producing offspring
#' @param fathers data frame of males producing offspring
#' @param clutch integer vector of length nrow(mothers) with clutch size per mother
#' @param probs numeric vector of length nrow(fathers) with probability of being a father
#' @param nstocks integer. number of stocks
#' @param Year integer. current year
cppFunction('DataFrame rcpp_mating(DataFrame mothers, DataFrame fathers, NumericVector clutch, 
  NumericVector probs, int nstocks, int Year) {
    
  int noffspring = sum(clutch);
  IntegerVector nos (nstocks+1); // for keeping track of number of offspring for each stock
  
  CharacterVector Dad (noffspring);
  CharacterVector Mum (noffspring);
  IntegerVector Stock (noffspring);
  
  // vectors for mother IDs and stock
  CharacterVector   MeM = mothers["Me"];
  IntegerVector   StockM = mothers["Stock"];
  
  // vectors for father IDs, stock and age
  CharacterVector   MeD = fathers["Me"];
  IntegerVector   StockD = fathers["Stock"];
  IntegerVector AgeD = fathers["AgeLast"];
  
  int iter = 0; int s; int m;
 
  for(s = 1; s <= nstocks; s++) {
  
    // subset mothers
    CharacterVector MeMS = MeM[StockM == s];
    if(MeMS.size() == 0) continue;

    // subset fathers
    CharacterVector   MeDS = MeD[StockD == s];
    IntegerVector  AgeDS = AgeD[StockD == s];
    if(MeDS.size() == 0) continue;
    NumericVector probsS = probs[StockD == s];
    
    // subset clutch
    NumericVector clutchInStock = clutch[StockM == s];
    int cs = sum(clutchInStock);
    nos[s] = cs;
    
    // assign fathers
    CharacterVector DadDS = sample(MeDS, cs, true, probsS);
    // now get the indexing right for the offsprings table
    for(m = 0; m < MeMS.size(); m++) {
      for(int i = iter; i < iter + clutchInStock[m]; i++){
        Mum[i] = MeMS[m];
        Dad[i] = DadDS[i-sum(nos[seq(0,s-1)])];
        Stock[i] = s;
      }
      iter = iter + clutchInStock[m];
    }
  }
  // fill in the rest
  CharacterVector Me(noffspring, NA_STRING);
  CharacterVector Sex(noffspring, "F");
  IntegerVector BirthY(noffspring, Year);
  IntegerVector DeathY(noffspring, NA_INTEGER);
  IntegerVector AgeLast(noffspring, 1);
  IntegerVector SampY(noffspring, NA_INTEGER);
  
  // if no mothers or fathers in a stock, it assigned stock 0
  // so omit those from result, if we have any
  DataFrame offspring = DataFrame::create(Named("Me") = Me[Stock != 0],
                                          Named("Sex") = Sex[Stock != 0],
                                          Named("Dad") = Dad[Stock != 0],
                                          Named("Mum") = Mum[Stock != 0],
                                          Named("BirthY") = BirthY[Stock != 0],
                                          Named("DeathY") = DeathY[Stock != 0],
                                          Named("Stock") = Stock[Stock != 0],
                                          Named("AgeLast") = AgeLast[Stock != 0],
                                          Named("SampY") = SampY[Stock != 0]);
  
  return offspring;
}')




# ---------------------------
# rcpp_mate_cmd

# rcpp version of Compound Multinomial Dirichlet sampling function
# called by cmate_cmd R function


#' Rcpp CMD function. The function loops over stocks and mother ages, then for each stock and mother age tabulates clutch size and then loops over unique clutch sizes to get N, a_dot and lambda for each clutch size, finally loops over mothers and samples fathers
#' this function is called within the R cmate_cmd function
#' the more potential fathers we have, the slower this thing is. breaking it up into stocks increases the speed significantly  
#' Note: rather than working with the UUIDs, we are working with integer vectors representing each individual in the mothers and fathers dataframe, then when we come back to R we match the integer vectors back to the individual's UUID because they are in the same order as they are in the mothers and fathers dfs
#' @param Nm number of mothers
#' @param Nm number of fathers
#' @param Nm integer vector of clutch size for each mother
#' @param probs numeric vector of relative reproductive output for each other (based on fecundity at age)
#' @param MS integer vector of stock ID for each mother
#' @param DS integer vector of stock ID for each father
#' @param Ma integer vector. mother ages
#' @param nstocks integer. max number of stocks
#' @param mating_sys string. one of "polyandry", "monogamy", or "mix"

cppFunction('DataFrame rcpp_mate_cmd(int Nm, int Nd, IntegerVector clutch, NumericVector probs, 
  NumericVector MS, NumericVector DS, IntegerVector Ma, int nstocks, std::string mating_sys) {
  
  int noffspring = sum(clutch);
  
  IntegerVector Dad (noffspring);
  IntegerVector Mum (noffspring);
  IntegerVector Stock (noffspring);
  
  IntegerVector fathers = seq(1,Nd);
  IntegerVector mothers = seq(1,Nm);
  
  double lambda;
  // parameters for scaling lambda by N and age, for monogamous mating
  double int_int = 2.7;
  double int_slope = -0.6;
  double slope_int = -0.4;
  double slope_slope = 0.2;
  
  int iter = 0;
  
  // loop over stocks
  for(int s = 1; s <= nstocks; s++){
  
    IntegerVector mothersS = mothers[MS == s]; // mothers in stock s
    if(mothersS.size() == 0) continue;
    IntegerVector fathersS = fathers[DS == s]; // fathers in stock s
    if(fathersS.size() == 0) continue;
    IntegerVector clutchS = clutch[MS == s]; // clutch sizes of mothers in stock s
    NumericVector probsS = probs[DS == s];
    IntegerVector MaS = Ma[MS == s]; // ages of mothers in stock
    int mmax_age = max(MaS); // max mother age
    
    NumericVector pp = probsS / sum(probsS); // probability of being a father
    
    // loop over mother ages
    for(int ma = 2; ma <= mmax_age; ma++){
       
       IntegerVector clutchSa = clutchS[MaS == ma]; // clutch sizes of mothers of age ma in stock s
       IntegerVector mothersSa = mothersS[MaS == ma]; // mothers in stock s of age ma
       if(mothersSa.size() == 0) continue;
       
       IntegerVector ct = as<IntegerVector>(table(clutchSa)); // count of clutchSa
       IntegerVector clutch_names = unique(clutchSa); // clutch sizes of count of clutchS
       clutch_names.sort(); // put them in increasing order
       
       int cstart; // first iteration of clutch size loop; this will depend on whether clutch size of 1 is present
      
      // clutch size 1 (if there is one)
      if(clutch_names[0] == 1){
        IntegerVector c1_dad = sample(fathersS, ct[0], true, probsS);
        IntegerVector idx = seq(iter,(iter + ct[0]-1));
        Dad[idx] = c1_dad;
        // no easy which implementation... need to get index of which clutchSa are 1
        IntegerVector tmpc = clutchSa[clutchSa == 1];
        IntegerVector clutchidx(tmpc.size()); 
        int iterc = 0;
        for(int i = 0; i < clutchSa.size(); i++){
          if(clutchSa[i] == 1) {
            clutchidx[iterc] = i;
            iterc ++;
          }
        }
        Mum[idx] = mothersSa[clutchidx];
        NumericVector tmpS1 (idx.size(), s);
        Stock[idx] = tmpS1;
        iter = iter + ct[0];
        cstart = 1;
      } else {
        cstart = 0;
      }
    
    
      // loop over clutch sizes. call them r
      for(int r = cstart; r < ct.size(); r++){ // start with clutch
        int N = clutch_names[r];
        double lam_int = int_int + N * int_slope;
        double lam_slope = slope_int + N * slope_slope;
        
        if (mating_sys == "monogamy") {
           //different lambda values depending on mother age
          if(ma == 2){
            lambda = 1.0001;
          } else {
            lambda = lam_int + ma * lam_slope;
            if (lambda > N-(N*0.03)) lambda = N-(N*0.03);
          }
        } else if (mating_sys == "polyandry"){
          lambda = 1.0001;
        } else { // mixed mating
          if(N == 2){
            lambda = 1.25;
          }
          else {
            lambda = N/2;
          }
        }
        
        double a_dot = (N - lambda) / (lambda - 1);
        NumericVector a = pp * a_dot;
      
        // another which. this time for which mothers have clutch size N
        IntegerVector tmpcr = clutchSa[clutchSa == N];
        IntegerVector clutchridx(tmpcr.size()); 
        int iterrc = 0;
        for(int i = 0; i < clutchSa.size(); i++){
          if(clutchSa[i] == N) {
            clutchridx[iterrc] = i;
            iterrc ++;
          }
        }
        IntegerVector clutchmothers = mothersSa[clutchridx];
        // loop over mothers
        for (int m = 0; m < ct[r]; m++){
            int tmpct = 0;
            NumericVector d(a.size());
            // every now and then, the d vector will be zero due to lambda
            // so here we give it 10 tries to get non-zero d vector
            // if by the 10th try we dont succeed the function will crash. just to make sure we dont get stuck in an infinite loop
            while(sum(d) == 0){  
                for(int j = 0; j < a.size(); j++){
                d[j] = rgamma(1, a[j], 1.0)[0];
              }
              d = d / sum(d);
              tmpct += 1;
              if (tmpct == 10)
                break;
            }
            int k = d.size();
            IntegerVector nmdraws(k);
            rmultinom(N, d.begin(), k, nmdraws.begin());
            IntegerVector samp = nmdraws[nmdraws > 0]; // number of offspring for each father
            // another which this time for which fathers are sampled
            IntegerVector sampled(N); 
            int iters = 0;
            for(int i = 0; i < nmdraws.size(); i++){
              if(nmdraws[i] > 0) {
                sampled[iters] = i; // index of fathers for this mothers offspring
                iters ++;
              }
            }
            IntegerVector didx = seq(iter, iter+N-1); // index of offspring that have current mother
            // now loop over sampled fathers for this mother m and repeat each one by the number of offspring he had with her 
            int tempiter = iter;
            for(int i = 0; i < samp.size(); i ++){
              IntegerVector tidx = seq(tempiter, tempiter + samp[i]-1);
              IntegerVector dadrep = rep(fathersS[sampled[i]], samp[i]);
              Dad[tidx] = dadrep;
              tempiter += samp[i];
            }
            IntegerVector mumrep = rep(clutchmothers[m], N);
            Mum[didx] = mumrep;
            NumericVector tmpS (didx.size(), s);
            Stock[didx] = tmpS; 
            iter = iter + N;
        } // m loop
      } // r loop
    } // mother ages loop
  } // stock loop
  
  DataFrame offspring = DataFrame::create(Named("Dad") = Dad[Stock != 0],
                                          Named("Mum") = Mum[Stock != 0],
                                          Named("Stock") = Stock[Stock != 0]);
  
  return offspring;
  
}')



# ---------------------------
# rcpp_mort

#' rcpp version of the mort2 function. faster
#' Note that nothing is returned from this function due to Rcpp's shallow copy property
#' @param indiv individuals data frame
#' @param year integer. current year
#' @param pdeath numeric vector. age-specific probability of dying
#' @param maxAge integer. plus group age
cppFunction('void rcpp_mort(DataFrame indiv, int year, NumericVector pdeath, int maxAge) {
  
  IntegerVector Age = indiv["AgeLast"];
  IntegerVector DeathY = indiv["DeathY"];
  
  for(int a = 1; a <= maxAge;a++) {
    IntegerVector AgeSub = Age[Age == a];
    int n = AgeSub.size();
    IntegerVector DeathYSub (n, NA_INTEGER);
    NumericVector probs = runif(n);

    DeathYSub[probs < pdeath[a-1]] = year;
    DeathY[Age == a] = DeathYSub;
  }
  
}')


# ---------------------------
# rcpp_mort_source

#' rcpp function to assign mortality source
#' @param indiv Individuals dataframe
#' @param probs numeric matrix of ages x mortality sources with relative mortality contribitions by age and source
#' @param dead binary integer vector of nrow(indiv) where 1s need mortality source assigned
#' @maxAge integer, maximum age of population
cppFunction('void rcpp_mort_source(DataFrame indiv, NumericMatrix probs, IntegerVector dead, int maxAge) {
  
  IntegerVector Age = indiv["AgeLast"];
  IntegerVector MortSource = indiv["MortSource"];
  IntegerVector AgeDead = Age[dead == 1];
  IntegerVector MortSourceDead = MortSource[dead == 1];
  IntegerVector Options = seq(1,3);
  
  for(int a = 1; a <= maxAge;a++) {
    IntegerVector AgeSub = AgeDead[AgeDead == a];
    if(AgeSub.size() == 0) continue;
    NumericVector probsAge = probs(a-1, _);
    int n = AgeSub.size();
    
    IntegerVector MortSourceSub = sample(Options, n, true, probsAge);
    MortSourceDead[AgeDead == a] = MortSourceSub;
  }
  MortSource[dead == 1] = MortSourceDead;
  
}')


# ---------------------------
# archive_dead

# not keeping track of age 1s that died from natural mortality to keep the df of dead
# individuals from getting too big
# since everyone has an equal probability of dying or getting sampled, this won't affect anything

archive_dead2 <- function (indiv = mort(), archive = make_archive()) 
{
  archive <- rbind(archive, indiv[!is.na(indiv[, 6]) & 
                                    !(indiv[, 8] == 1 & (is.na(indiv[,11]) | indiv[,11] == 3)), ])
  return(archive)
}



# ---------------------------
# birthdays

# modified the birthdays function to make the final age a plus group

birthdays2 <- function (indiv = makeFounders()) 
{
  indiv[is.na(indiv[, 6]) & indiv[, 8] < maxAge, 8] <- indiv[is.na(indiv[, 6]) & indiv[, 8] < maxAge, 8] + 
    1L
  return(indiv)
}


# ---------------------------
# fem_stock_switch

#' function for switching stock
#' @param indiv individuals dataframe
#' @param switch_prob number between 0 and 1 spawning site fidelity

fem_stock_switch <- function(indiv = makeFounders(), switch_prob){
  # distance matrix for switching stocks
  switch_mat <- 1/as.matrix(dist(seq(from = 1, by = 8, length.out = nstocks), upper = T))
  diag(switch_mat) <- 0
  switch_mat <- sweep(switch_mat, 1, FUN = "/", rowSums(switch_mat))
  # males stay the stock they are, females may switch
  fa <- which(indiv$Sex == "M") # index males - they don't switch
  sw <- rbinom(nrow(indiv), 1, (1 - switch_prob)) # binomial probability of switching stock
  sw[fa] <- 0 # males not switching stocks
  cs <- indiv$Stock[sw == 1] # current stock of individuals switching
  ns <- vector("numeric", length = length(cs))
  si <- 1:nstocks
  if (nstocks > 2){
    for (i in si) {
      ind <- which(cs == i)
      #ns[ind] <- sample(si[-which(si == i)], length(ind), replace = TRUE)
      ns[ind] <- sample(si, length(ind), replace = TRUE, prob = switch_mat[i,])
    }
    indiv[sw == 1, 7] <- ns
  } else {
    indiv[sw == 1, 7] <- ifelse(indiv[sw == 1, 7] == 1, 2, 1) # perform stock switch
  }
  return(indiv)
}



