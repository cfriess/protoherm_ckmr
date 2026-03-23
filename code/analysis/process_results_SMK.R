

library(tidyverse)
library(nplyr)

# snakemake harvested outputs file for keyObjs:
res_raw <- readRDS("file_path/results.rds")

colnames(res_raw)[colnames(res_raw) == "rds"] <- "data"
res_raw$mating[res_raw$mating == "cmd"] <- "monogamy"
res_raw$mating[res_raw$mating == "random"] <- "polyandry"
res_raw$id <- 1:nrow(res_raw)
res_raw$sample_by_sex <- as.numeric(ifelse(res_raw$sample_by_sex == "TRUE", "1", res_raw$sample_by_sex))


res <- res_raw %>%
  mutate(scenario = case_when(
    sample_by_sex %in% 2:3 ~ "Suppl Scen Male & Female Sampling",
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 1 & 
      hzt == 1 & hze == 1 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3 ~ "Base", 
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 1 & 
      hzt == 0.9 & hze == 1 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3 ~ "sHmt < eHmt",
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 1 & 
      hzt == 1 & hze == 0.9 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3 ~ "eHmt < sHmt",
    mating == "polyandry" & eofsr == 3 & ssf == 0.01 & fec_known == 1 & 
      hzt == 1 & hze == 1 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3 ~ "skipped spawning",
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 2 & 
      hzt == 1 & hze == 1 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3  ~ "mRV steeper",
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 3 & 
      hzt == 1 & hze == 1 & nll_type == 1 & et %in% 0:3 &
      ckmr_ssmult == 1.0 & nsampyrs == 3  ~ "mRV flatter",
    mating == "polyandry" & eofsr == 0 & ssf == 0.01 & fec_known == 1 & 
      hzt == 1 & hze == 1 & nll_type == 1 & et %in% 0:3 ~ "Suppl Sampling Intensity",
    TRUE ~ "not used"
  )) %>%
  filter(scenario != "not used") %>%
  mutate(scen_type = ifelse(scenario %in% c("Suppl Sampling Intensity",
                                            "Suppl Scen Male & Female Sampling"), "Supplement", "Main"))
table(res$scenario[res$scen_type == "Main"]) # 1600 per scenario

# N for sample size and # of sampling year runs

res <- res %>%
  mutate(
    et = case_when(
      et == 0 ~ "fixed trans & mRV",
      et == 1 ~ "est trans & fixed mRV",
      et == 2 ~ "fixed trans & est mRV",
      et == 3 ~ "est trans & mRV"
    ),
    fec_known = case_when(
      fec_known == 1 ~ "mRV = fRV",
      fec_known == 2 ~ "mRV steeper",
      fec_known == 3 ~ "mRV flatter"
    ),
     hze = ifelse(hze == 0.9, "hze = 0.9", "hze = 1"),
     hzt = ifelse(hzt == 0.9, "hzt = 0.9", "hzt = 1"),
     eofsr = case_when(
       eofsr == 0 ~ "No Sperm Lim",
       eofsr == 1 ~ "Prop. Red. Fert",
       eofsr == 2 ~ "Prop. Skipped Spawning",
       eofsr == 3 ~ "Disprop. Skipped Spawning"
     ),
     ssf = ifelse(ssf == 0.01, "No Fem SS Fidel", "80% Fem SS Fidel"),
     nsampyrs = ifelse(nsampyrs == 3, "03 Yrs", "10 Yrs"),
     ckmr_ssmult = case_when(
       ckmr_ssmult == 0.5 ~ "50%",
       ckmr_ssmult == 1 ~ "100%",
       ckmr_ssmult == 1.5 ~ "150%"
     ),
     nll_type = ifelse(nll_type == 0, "HSP", "HSP + POP"),
    sample_by_sex = case_when(
      sample_by_sex == 1 ~ "age 1+ sex ratio",
      sample_by_sex == 2 ~ "age 2+ sex ratio",
      sample_by_sex == 3 ~ "females only sampled"
    )
  ) %>%
  mutate(om_sc = case_when(
    eofsr %in% "Disprop. Skipped Spawning"  ~ "skipped spawning",
    fec_known == "mRV flatter" ~ "mRV flatter",
    fec_known == "mRV steeper" ~ "mRV steeper",
    fec_known == "mRV = fRV" ~ "Base",
    TRUE ~ "check me"
  ),
  em_sc = case_when(
    hzt == "hzt = 0.9" ~ "sim Hmt < est Hmt",
    hze == "hze = 0.9" ~ "est Hmt < sim Hmt",
    nll_type == "HSP" ~ "HSPs only",
    TRUE ~ "Base"
  ))

# N for the main paper runs:
(N <- nrow(res %>% filter(scen_type == "Main"))) # 9600
nrow(res %>% filter(scen_type == "Supplement")) # 54,400

res_compn <- res %>% 
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, data, id, scenario, om_sc, em_sc, scen_type) %>%
  mutate(abund_res = map(data, "abund_res"),
         admb_cormat = map(data, "admb_cor"),
         para_res = map(data, "para_res"),
         sim_pop = map(data, "sim_pop"),
         ckmr_obs = map(data, "ckmr_obs"),
         ckmr_true = map(data, "ckmr_true"),
         admb_exit = map_dbl(data, "admb_exit"),
         max_grad = map_dbl(data, "max_grad"),
         admb_run_res = map(data, "admb_run_res")
         ) 

### admb failures
admb_fail_ids <- unique(res_compn$id[res_compn$admb_exit == 1])

### boundary estimation

rbar_bounds <- exp(c(6,17))
f_bounds <- exp(c(-4.6,0.5))
tsd_bounds <- exp(c(0.9,2.4))
mexp_bounds <- exp(c(-4,2))

paras <- res_compn %>% 
  filter(admb_exit == 0) %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, para_res, om_sc, em_sc, scen_type) %>%
  unnest(para_res) %>%
  mutate(lb = case_when(
    par_name == "log_rbar" ~ rbar_bounds[1],
    par_name == "log_tsd" ~ tsd_bounds[1],
    par_name == "log_mfec_exp" ~ mexp_bounds[1],
    par_name %in% c("log_fbarc", "log_fbarr") ~ f_bounds[1]
  ),
  ub = case_when(
    par_name == "log_rbar" ~ rbar_bounds[2],
    par_name == "log_tsd" ~ tsd_bounds[2],
    par_name == "log_mfec_exp" ~ mexp_bounds[2],
    par_name %in% c("log_fbarc", "log_fbarr") ~ f_bounds[2]
  )) %>%
  mutate(bounds_flag = ifelse(value <= lb + 0.01 * lb | value > ub - 0.01 * ub, 1, 0))
table(paras$bounds_flag, paras$par_name) # some bounds hit for male RV and some for tsd, some for log_rbar

length(unique(paras$id))
table(paras$par_name)

bounds_exclude_ids <- unique(paras$id[paras$bounds_flag == 1])

### high max gradients

max_grad_thresh <- 0.1 # set to 1e8 for no gradient filtering or 0.1 for more restrictive gradient filtering

high_grad_ids <- unique(res_compn$id[res_compn$max_grad > max_grad_thresh & !is.na(res_compn$max_grad)])
high_grad_ids <- high_grad_ids[!high_grad_ids %in% bounds_exclude_ids]

exclude_ids <- unique(c(admb_fail_ids, bounds_exclude_ids, high_grad_ids))

res_compn$bounds_fail <- 0
res_compn$bounds_fail[res_compn$id %in% bounds_exclude_ids] <- 1

res_compn$gradient_fail <- 0
res_compn$gradient_fail[res_compn$id %in% high_grad_ids] <- 1

table(res_compn$admb_exit, res_compn$scen_type) # 227 out of 9600 failed
table(res_compn$bounds_fail, res_compn$scen_type) # 779 out of 9373 failed
table(res_compn$gradient_fail, res_compn$scen_type) # 38 out of 8821 failed

res_compn$exclude <- 0
res_compn$exclude[res_compn$id %in% exclude_ids] <- 1

table(res_compn$exclude, res_compn$scen_type) # 8556 out of 9600 used

sim_pars <- paras %>%
  filter(!id %in% exclude_ids) %>%
  mutate(perror = (value - sim_value)/sim_value * 100,
         aperror = abs(value - sim_value)/sim_value * 100) %>%
  bind_rows(paras %>% filter(par_name == "log_tsd") %>%
              mutate(et = "fixed trans & mRV", 
                     perror = 0,
                     aperror = 0,
                     value = sim_value)) %>%
  bind_rows(paras %>% filter(par_name == "log_tsd") %>%
              mutate(et = "fixed trans & est mRV", 
                     perror = 0,
                     aperror = 0,
                     value = sim_value)) %>%
  bind_rows(paras %>% filter(par_name == "log_mfec_exp") %>%
              mutate(et = "est trans & fixed mRV", 
                     perror = 0,
                     aperror = 0,
                     value = 1.4935)) %>%
  bind_rows(paras %>% filter(par_name == "log_mfec_exp") %>%
              mutate(et = "fixed trans & mRV", 
                     perror = 0,
                     aperror = 0,
                     value = 1.4935)) %>%
  mutate(perror = (value - sim_value)/sim_value * 100,
         aperror = abs(value - sim_value)/sim_value * 100) %>%
  mutate(par_name = case_when(
    par_name == "log_rbar" ~ "R", 
    par_name == "log_fbarc" ~ "Fcom",  
    par_name == "log_fbarr" ~ "Frec",  
    par_name == "log_tsd" ~ "TFsd",
    par_name == "log_mfec_exp" ~ "mRV_exp"
  ))


saveRDS(sim_pars, "data/results/full_sim_summary_results/sim_pars_final.rds", compress = "xz")

admb_res <- res_compn %>%
  filter(!id %in% exclude_ids) %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, admb_run_res, om_sc, em_sc, scen_type) %>%
  unnest(admb_run_res)
run_cnt <- admb_res %>% 
  filter(!is.na(admb_run_cnt)) %>%
  group_by(id, scen_type) %>%
  summarize(max_runs = max(admb_run_cnt),
            best_nll = which.min(nll)) 
table(run_cnt$max_runs, run_cnt$best_nll, run_cnt$scen_type)


abund_res <- res_compn %>%
  filter(!id %in% exclude_ids) %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, abund_res, om_sc, em_sc, scen_type) %>%
  unnest(abund_res) %>%
  mutate(perror = (est_value - sim_value)/sim_value * 100,
         aperror = abs(est_value - sim_value)/sim_value * 100) %>%
  filter(data!= "Recruits") 

# this is a big file, use maximal compression for github, but read and write speed will be lower
saveRDS(abund_res, "data/results/full_sim_summary_results/abund_res_final.rds", compress = "xz")

sim_pop <- res_compn %>%
  #filter(admb_exit == 0) %>%
  filter(ckmr_seed == 1, et == "fixed trans & mRV", nll_type == "HSP + POP", 
         hzt == "hzt = 1", hze == "hze = 1", sample_by_sex == "age 1+ sex ratio",
         nsampyrs == "03 Yrs", ckmr_ssmult == "100%") %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known,
                id, scenario, sim_pop) %>%
  unnest(sim_pop)

saveRDS(sim_pop, "data/results/full_sim_summary_results/sim_pop_final.rds")

ckmr_samps <- res_compn %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, ckmr_true, exclude,
                om_sc, em_sc, scen_type) %>%
  nest_summarize(ckmr_true, 
                 POD = sum(PO),
                 POS = sum(MO),
                 HSD = sum(HSD),
                 HSS = sum(HSS),
                 Not = sum(Not)) %>%
  unnest(ckmr_true) 

saveRDS(ckmr_samps, "data/results/full_sim_summary_results/ckmr_samps_final.rds")

ckmr_obs_samps <- res_compn %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, ckmr_obs, exclude,
                om_sc, em_sc, scen_type) %>%
  nest_summarize(ckmr_obs, 
                 POD = sum(PO),
                 POS = sum(MO),
                 HSD = sum(HSD),
                 HSS = sum(HSS),
                 Not = sum(Not)) %>%
  unnest(ckmr_obs)

saveRDS(ckmr_obs_samps, "data/results/full_sim_summary_results/ckmr_obs_samps_final.rds")


samp_summary <- res_compn %>%
  mutate(samp_summary = map(data, "ckmr_samps")) %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, samp_summary, exclude,
                om_sc, em_sc, scen_type) %>%
  nest_mutate(samp_summary, LifeStage = ifelse(AgeLast == 1, "Juv", "Adult")) %>%
  nest_group_by(samp_summary, Sex, LifeStage) %>%
  nest_summarize(samp_summary, N = sum(N)) %>%
  nest_ungroup(samp_summary) %>%
  unnest(samp_summary)

saveRDS(samp_summary, "data/results/full_sim_summary_results/ckmr_samp_summary_final.rds")




#### Parameter corrleation

cor_res_et <- res_compn %>%
  filter(!id %in% exclude_ids, et == "est trans & mRV") %>%
  dplyr::select(ped_rep, eofsr, mating, ssf, fec_known, sample_by_sex,
                nsampyrs, ckmr_ssmult, et, nll_type, hzt, hze,
                ckmr_seed, id, scenario, admb_cormat, om_sc, em_sc,
                scen_type) 

cor_res_et$rbar_tsd_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
  mat[1, 4]
})

cor_res_et$rbar_exp_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
  mat[1, 5]
})

cor_res_et$tsd_exp_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
  mat[4, 5]
})

# cor_res_et$rbar_Fcom_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
#   mat[1, 2]
# })
# cor_res_et$rbar_Frec_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
#   mat[1, 3]
# })
# cor_res_et$Fcom_Frec_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
#   mat[2, 3]
# })
# cor_res_et$Fcom_tsd_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
#   mat[2, 4]
# })
# cor_res_et$Frec_tsd_cor <- sapply(cor_res_et$admb_cormat, function(mat) {
#   mat[3, 4]
# })
# 
#
cor_res_et$ckmr_ssmult <- factor(cor_res_et$ckmr_ssmult, levels = c("50%","100%", "150%"))


saveRDS(cor_res_et, "data/results/full_sim_summary_results/cor_res_et_final.rds")


