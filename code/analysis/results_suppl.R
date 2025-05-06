
## ---------------------------
##
## Script name: results_suppl.R
##
## Purpose of script: generate figures for study results
##
## Author: Claudia & Eric
##
## Date Created: 2025-03-20
##
## Notes:
##
##
##
## ---------------------------

library(tidyverse)
library(nplyr)
library(patchwork)
library(ggh4x)
library(cowplot)
library(gridExtra)
library(grid)
library(scales)
library(latex2exp)
library(patchwork)
library(here)

theme_set(theme_minimal_hgrid())



# ---------------------------------------------
#
# data prep
#
# ---------------------------------------------

path <- "data/results/full_sim_summary_results"

source(here("code", "main_sim", "sim_inputs.R")) # simulation inputs


### ----------------------------------
# read in results files

et_levels <- c("fixed trans & mRV", 
               "fixed trans & est mRV", "est trans & fixed mRV", "est trans & mRV")
om_sc_levels <- c("Base", "skipped spawning", "mRV flatter", "mRV steeper")
em_sc_levels <- c("Base", "HSPs only", "est Hmt < sim Hmt", "sim Hmt < est Hmt")

# parameter correlations
cor_res_et <- readRDS(file.path(here(), path, "cor_res_et_final.rds"))
cor_res_et$om_sc <- factor(cor_res_et$om_sc, levels = om_sc_levels)

# ckmr sampling summary
samp_summary <- readRDS(file.path(here(), path, "ckmr_samp_summary_final.rds")) 
samp_summary$om_sc <- factor(samp_summary$om_sc, levels = om_sc_levels)
samp_summary$ckmr_ssmult <- factor(samp_summary$ckmr_ssmult, levels = c("50%","100%", "150%"))
samp_summary$et <- factor(samp_summary$et, levels = et_levels)

# abundance estimate results
abund_res <- readRDS(file.path(here(), path, "abund_res_final.rds")) %>%
  mutate(squared_diff = (sim_value - est_value)^2)
abund_res$ckmr_ssmult <- factor(abund_res$ckmr_ssmult, levels = c("50%","100%", "150%"))
abund_res$et <- factor(abund_res$et, levels = et_levels)
abund_res$om_sc <- factor(abund_res$om_sc, levels = om_sc_levels)
abund_res$em_sc <- factor(abund_res$em_sc, levels = em_sc_levels)

# simulated parameters
sim_pars <- readRDS(file.path(here(), path, "sim_pars_final.rds"))
sim_pars$ckmr_ssmult <- factor(sim_pars$ckmr_ssmult, levels = c("50%","100%", "150%"))
sim_pars$et <- factor(sim_pars$et, levels = et_levels)
sim_pars$om_sc <- factor(sim_pars$om_sc, levels = om_sc_levels)
sim_pars$em_sc <- factor(sim_pars$em_sc, levels = em_sc_levels)
sim_pars$scenario <- as.factor(sim_pars$scenario)
sim_pars$scenario <- factor(sim_pars$scenario, levels(sim_pars$scenario)[c(22,18:21,1,10:17,2:9)])

sim_pars <- sim_pars %>%
  mutate(par_name_long = case_when(
    par_name == "R" ~ "Mean Recruitment",
    par_name == "Fcom" ~ "Mean Commercial F",  
    par_name == "Frec" ~ "Mean Recreational F",  
    par_name == "TFsd" ~ "Sex Transition Fct SD",
    par_name == "mRV_exp" ~ "Male RV Fct Exp"
  ))


# true ckmr samples
ckmr_samps <- readRDS(file.path(here(), path, "ckmr_samps_final.rds")) 
ckmr_samps$ckmr_ssmult <- factor(ckmr_samps$ckmr_ssmult, levels = c("50%","100%", "150%"))
ckmr_samps$et <- factor(ckmr_samps$et, levels = et_levels)
ckmr_samps$om_sc <- factor(ckmr_samps$om_sc, levels = om_sc_levels)
ckmr_samps$em_sc <- factor(ckmr_samps$em_sc, levels = em_sc_levels)

ckmr_samps_long <- ckmr_samps %>%
  gather(pair_type, count, POD:Not) %>%
  mutate(pair_type = ifelse(pair_type == "Not", "Unrelated", pair_type))
ckmr_samps_long$pair_type <- factor(ckmr_samps_long$pair_type, levels = c("POD", "POS", "HSD", "HSS", "Unrelated"))

# observed ckmr samples
ckmr_obs_samps <- readRDS(file.path(here(), path, "ckmr_obs_samps_final.rds"))
ckmr_obs_samps$ckmr_ssmult <- factor(ckmr_obs_samps$ckmr_ssmult, levels = c("50%","100%", "150%"))
ckmr_obs_samps$et <- factor(ckmr_obs_samps$et, levels = et_levels)
ckmr_obs_samps$om_sc <- factor(ckmr_obs_samps$om_sc, levels = om_sc_levels)
ckmr_obs_samps$em_sc <- factor(ckmr_obs_samps$em_sc, levels = em_sc_levels)

ckmr_obs_samps_long <- ckmr_obs_samps %>%
  gather(pair_type, count, POD:Not) %>%
  mutate(pair_type = ifelse(pair_type == "Not", "Unrelated", pair_type))
ckmr_obs_samps_long$pair_type <- factor(ckmr_obs_samps_long$pair_type, levels = c("POD", "POS", "HSD", "HSS", "Unrelated"))

# base case levels
base_eofsr <- "No Sperm Lim"
base_ssf <- "No Fem SS Fidel"
base_mating <- "polyandry"
base_ckmr_nsampyrs <- "03 Yrs"
base_ckmr_ssmult <- "100%"
base_hze <- "hze = 1"
base_hzt <- "hzt = 1"
base_nll <- "HSP + POP"
base_fec <- "Base"
base_et <- "fixed trans & mRV"
samp_sex_base <- "age 1+ sex ratio"


# -----------------------------
# define some common quantities

get_perror = function(sim, est) (est - sim)/sim * 100

ytitle <- ggdraw() + draw_label("Percent Relative Error", 
                                size = 16, 
                                fontface = "plain", 
                                angle = 90)

xtitle <- ggdraw() + draw_label("CKMR Sample Size", 
                                size = 16, 
                                fontface = "plain")

fyrplot <- 75
prop2plot <- 0.2

em_labels <- c(
  "fixed trans & mRV" = expression("fixed" ~ italic(P)^{"F" %->% "M"} ~ " & RV"["M"]),
  "fixed trans & est mRV" = expression("fixed" ~ italic(P)^{"F" %->% "M"} ~ ", free RV"["M"]),
  "est trans & fixed mRV" = expression("free" ~ italic(P)^{"F" %->% "M"} ~ ", fixed RV"["M"]),
  "est trans & mRV" = expression("free" ~ italic(P)^{"F" %->% "M"} ~ " & RV"["M"])
)

om_labels <- c(
  "Base" = "Base",
  "skipped spawning" = "Skipped spawning",
  "mRV flatter" = expression(RV[M]^{"OM"} ~ "flatter than" ~ RV[F]^{"OM"}),
  "mRV steeper" = expression(RV[M]^{"OM"} ~ "steeper than" ~ RV[F]^{"OM"})
)




###############################
#
# Supplemental Figures
#
###############################

plot_age <- seq(1,33,0.1)
plot_age_Male <- seq(4,33,0.1)

#--------------------------------
# Figure S1: age schedules

# maturity at age
female_mat = exp(m1+m2*plot_age)/(1+exp(m1+m2*plot_age))
female_mat[plot_age <= 1] <- 0 # age 1 not mature

# probability of transitioning at age
trans_prob <- pnorm(q=plot_age, mean=t1, sd=t2)
trans_prob[plot_age <= 3] <- 0 # zero percent chance of transitioning for ages 1 to 3

Mage[1] <- NA
Mdf <- data.frame(Age = 1:maxAge, M = Mage)

age_schedule_df <- data.frame(Age = plot_age, FemaleMaturity = female_mat,
                              TransitionProb = trans_prob) %>%
  gather(Metric, Probability, -Age)

age_schedule_plot <- ggplot(age_schedule_df) +
  geom_line(aes(x = Age, y = Probability, linetype = Metric)) +
  geom_line(data = Mdf, aes(x = Age, y = M), linetype = "dashed") +
  theme_bw() + ylab("Probability/ Mortality") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background=element_blank(),
        legend.title = element_blank())


#--------------------------------
# Figure S2: Fishing mortality and vulnerability


### vulnerability

len_plot <- linf * (1-exp(-k*(plot_age-t0)))

vul_plot = (len_plot/lmax)^(lmax/p) * exp((lmax-len_plot)/p)
vul_plot[len_plot < 200] <- 0 # for now, no fishing on age 0
vul_plot <- vul_plot / max(vul_plot)

llvul_plot <- 1 / (1+exp(-a1*(plot_age-a0)))

vul_df <- data.frame(Age = plot_age, Recreational = vul_plot, Commercial = llvul_plot) %>%
  gather(Metric, Vulnerability, -Age)

vul_plot <- ggplot(vul_df) +
  geom_line(aes(x = Age, y = Vulnerability, linetype = Metric)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background=element_blank(),
        legend.title = element_blank())

#### fishing mortality
fplot_dat <- Fsim %>% gather(Fleet, F, -year)

fplot <- ggplot(fplot_dat) +
  geom_line(aes(x = year, y = F, linetype = Fleet)) +
  theme_bw() + ylab("Fishing Mortality") + xlab("Year") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.background=element_blank(),
        legend.title = element_blank())

F_vul_plot <- vul_plot + fplot +
  plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A')


#--------------------------------
# Figure S3: Number of individuals sampled for CKMR

nsamps_plot <- ggplot(samp_summary %>%
                        filter(LifeStage != "Juv") %>%
                        mutate(Sex = ifelse(Sex == "F", "Female", "Male")) %>%
                        filter(om_sc == "Base", em_sc == "Base", et == base_et,
                               sample_by_sex == samp_sex_base)) +
  geom_boxplot(aes(x = ckmr_ssmult, y = N, fill = nsampyrs)) +
  facet_wrap(~Sex, scales = "free_y") +
  xlab("CKMR Sample Size") + ylab("Number of Individuals Sampled") +
  guides(color = guide_legend(title = "Number of CKMR Sampling Years")) +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


#--------------------------------
# Figure S4: Number of observed kin pairs by CKMR sampling

ckmr_pairs_by_sampling <- ckmr_samps_long %>% filter(scenario == "S0")

ckmr_pairs_by_sampling_TEX <- ckmr_pairs_by_sampling
levels(ckmr_pairs_by_sampling_TEX$pair_type) <- c(
  POD = TeX("POP$^\\neq$"),
  POS = TeX("POP$^=$"),
  HSD = TeX("HSP$^\\neq$"),
  HSS = TeX("HSP$^=$"),
  Unrelated = TeX("Unrelated")
)

suppl_ckmr_pairs_plot <- ggplot(ckmr_pairs_by_sampling_TEX) +
  geom_boxplot(aes(x = ckmr_ssmult, y = count, fill = nsampyrs)) +
  facet_wrap(~pair_type, scales = "free", ncol = 5, labeller = label_parsed) +
  ylab("Count") +
  guides(color = guide_legend(title = "Number of CKMR Sampling Years")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size = 16)
  ) +
  xlab("CKMR Sample Size")

# ckmr_pairs_by_sampling %>% 
#   group_by(pair_type, ckmr_ssmult) %>%
#   summarize(medN = median(count)) %>%
#   spread(pair_type, medN)




#--------------------------------
# Figure S5: PRE for POPs + HSPs vs HSPs only for base and skipped spawning


err_hsps <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, om_sc == "skipped spawning") %>%
  group_by(om_sc, em_sc, et, data, id) %>%
  summarize(perror = median(perror))

em_hsps_p1 <- 
  ggplot(err_hsps %>% filter(em_sc %in%  "Base"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  ggtitle("Base EM (HSPs + POPs)") +
  ylab("Percent Relative Error") +
  coord_cartesian(ylim = c(-55,100)) +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

em_hsps_p2 <- 
  ggplot(err_hsps %>% filter(em_sc %in% "HSPs only"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  ggtitle("HSPs only") +
  ylab("Percent Relative Error") +
  coord_cartesian(ylim = c(-55,100)) +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "B")


ckmr_hsp_plot <- em_hsps_p1 + em_hsps_p2 +
  plot_layout(axis_titles = "collect",
              guides = "collect") &
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank()) &
  guides(fill = guide_legend(title = "EM Scenario", nrow = 2,
                             byrow = T))



#--------------------------------
# Figure S6: ckmr pairs by EM scenario


ckmr_pairs_by_scenario2 <- ckmr_obs_samps_long %>% 
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, ckmr_ssmult == base_ckmr_ssmult,
         om_sc == "Base") 

ckmr_pairs_by_scenario2_TEX <- ckmr_pairs_by_scenario2
levels(ckmr_pairs_by_scenario2_TEX$pair_type) <- c(
  POD = TeX("POP$^\\neq$"),
  POS = TeX("POP$^=$"),
  HSD = TeX("HSP$^\\neq$"),
  HSS = TeX("HSP$^=$"),
  Unrelated = TeX("Unrelated")
)

suppl_ckmr_pairs_plot2 <- ggplot(ckmr_pairs_by_scenario2_TEX) +
  geom_boxplot(aes(x = em_sc, y = count, fill = em_sc)) +
  facet_wrap(~pair_type, labeller = label_parsed,
             scales = "free", ncol = 3) +
  scale_fill_discrete(
    labels = c(
      "Base" = "Base",
      "HSPs only" = "HSPs only",
      "est Hmt < sim Hmt" = expression(H[mt]^{"EM"} < H[mt]^{"OM"}),
      "sim Hmt < est Hmt" = expression(H[mt]^{"OM"} < H[mt]^{"EM"})
    )
  ) +
  ylab("Number of Pairs") + xlab("") +
  theme(axis.text.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.75,0.25),
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        strip.text.x = element_text(size = 16)) +
  guides(fill = guide_legend(title = "EM Scenario"))



#--------------------------------
# Figure S7: estimated parameters by EM scenario

layout3 <- "
BCF
DEA
" 

pars_nll <- sim_pars %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, om_sc == "Base") 

levels(pars_nll$et)[levels(pars_nll$et) == "fixed trans & est mRV"] <- "fixed trans, free mRV"
levels(pars_nll$et)[levels(pars_nll$et) == "est trans & fixed mRV"] <- "free trans, fixed mRV"
levels(pars_nll$et)[levels(pars_nll$et) == "est trans & mRV"] <- "free trans & mRV"

par_est_nll <- pars_nll %>%
  split(., .$par_name_long) %>%
  map(function(x){
    p <-  ggplot(x, aes(x = et, y = perror)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
      geom_boxplot(aes(fill = em_sc),
                   outlier.size = 0.8, outlier.fill = NULL, 
                   outlier.shape = 21, outlier.stroke = 0.3) +
      scale_fill_discrete(
        labels = c(
          "Base" = "Base",
          "HSPs only" = "HSPs only",
          "est Hmt < sim Hmt" = expression(H[mt]^{"EM"} < H[mt]^{"OM"}),
          "sim Hmt < est Hmt" = expression(H[mt]^{"OM"} < H[mt]^{"EM"})
        )
      ) +
      facet_nested(. ~ par_name_long) +
      scale_x_discrete(labels = scales::wrap_format(15)) +
      ylab("Percent Relative Error") +
      guides(fill = guide_legend(title = "EM Scenario")) 
    
    if (unique(x$par_name) == "mRV_exp") {
      p <- p + coord_cartesian(ylim = c(-100, 200)) # 1 outliers removed incorp. 5% expand
    }
    
    if (unique(x$par_name) == "R") {
      p <- p + coord_cartesian(ylim = c(-73, 500)) # 1 outliers removed incorp. 5% expand
    }
    
    return(p)
    
  })

para_res_plot_nll <- wrap_plots(par_est_nll, design = layout3) + guide_area() + 
  plot_layout(axis_titles = "collect", guides = "collect", axes = "collect") & 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))




#--------------------------------
# Figure S8: estimated parameters by OM scenario

pars_mRV <- sim_pars %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, em_sc == "Base"
  ) 

levels(pars_mRV$et)[levels(pars_mRV$et) == "fixed trans & est mRV"] <- "fixed trans, free mRV"
levels(pars_mRV$et)[levels(pars_mRV$et) == "est trans & fixed mRV"] <- "free trans, fixed mRV"
levels(pars_mRV$et)[levels(pars_mRV$et) == "est trans & mRV"] <- "free trans & mRV"

par_est <- pars_mRV %>%
  split(., .$par_name_long) %>%
  map(function(x){
    p <-  ggplot(x, aes(x = et, y = perror)) +
      geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
      geom_boxplot(aes(fill = om_sc),
                   outlier.size = 0.8, outlier.fill = NULL, 
                   outlier.shape = 21, outlier.stroke = 0.3) +
      scale_fill_discrete(
        labels = om_labels
      ) +
      facet_nested(. ~ par_name_long) +
      scale_x_discrete(labels = scales::wrap_format(15)) +
      ylab("Percent Relative Error") +
      guides(fill = guide_legend(title = "OM Scenario")) 
    
    if (unique(x$par_name) == "mRV_exp") {
      p <- p + coord_cartesian(ylim = c(-90, 1400)) # 1 outliers removed incorp. 5% expand
    }
    
    if (unique(x$par_name) == "R") {
      p <- p + coord_cartesian(ylim = c(-72, 150)) # 1 outliers removed incorp. 5% expand
    }
    
    return(p)
    
  })

para_res_plot <- wrap_plots(par_est, design = layout3) + guide_area() + 
  plot_layout(axis_titles = "collect", guides = "collect", axes = "collect") & 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



#--------------------------------
# Figure S9: parameter correlations

cor_dat <- cor_res_et %>%
  filter(scenario %in% c(paste0("EM_S", c(5,7,9,17))), nsampyrs == base_ckmr_nsampyrs) %>%
  dplyr::select(nsampyrs, ckmr_ssmult, rbar_tsd_cor:tsd_exp_cor,fec_known, om_sc) %>%
  gather(comp, cor, rbar_tsd_cor:tsd_exp_cor) %>%
  mutate(comp = case_when(
    comp == "rbar_tsd_cor" ~ "Recruitment & Transition Fct",
    comp == "rbar_exp_cor" ~ "Recruitment & Male RV",
    comp == "tsd_exp_cor" ~ "Transition Fct & Male RV"
  ))


cor_plot <- ggplot(cor_dat) +
  geom_boxplot(aes(x = ckmr_ssmult, y = cor, fill = om_sc)) +
  facet_wrap(~comp) +
  scale_fill_discrete(
    labels = om_labels
  ) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  xlab("CKMR Sample Size") +
  guides(fill = guide_legend(title = "OM Scenario", nrow = 2, byrow = T)) +
  ylab("Correlation")



#--------------------------------
# Figures A10-A13: PRE and NRMSE comparison

error_metrics_scens <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult,
         !(om_sc == "skipped spawning" & em_sc == "HSPs only")) %>%
  mutate(sc3 = case_when(
    em_sc != "Base" ~ em_sc,
    om_sc == "Base" ~ "Base",
    TRUE ~ om_sc
  )) %>%
  bind_rows(err_ckmr_sampling <- abund_res %>%
              filter(scenario %in% c("S0","EM_S4", "EM_S5"),
                     !(ckmr_ssmult == base_ckmr_ssmult &
                         nsampyrs== base_ckmr_nsampyrs)) %>%
              mutate(sc3 = case_when(
                nsampyrs == "03 Yrs" & ckmr_ssmult == "50%" ~ "3 Yrs, 50% BaseN",
                nsampyrs == "03 Yrs" & ckmr_ssmult == "150%" ~ "3 Yrs, 150% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "50%" ~ "10 Yrs, 50% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "100%" ~ "10 Yrs, 100% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "150%" ~ "10 Yrs, 150% BaseN",
                TRUE ~ "check me"
              )))

error_metrics <- error_metrics_scens %>%
  group_by(sc3, et, data, id, ped_rep, ckmr_seed) %>% # summarize by id
  summarize(pre = median(perror),
            ape = median(aperror),
            sd_ape = sd(aperror),
            med_sim = median(sim_value),
            med_est = median(est_value),
            squared_diff = (med_est-med_sim)^2) %>%
  ungroup() %>%
  group_by(sc3, et, data) %>% # summarize by scenario
  summarize(med_PRE = median(pre),
            med_APE = median(ape),
            iqr_PRE = IQR(pre),
            N = length(unique(id)),
            meanobs = mean(med_sim),
            rmse = sqrt((1/N) * sum(squared_diff)),
            nrmse = rmse/meanobs,
            min_PRE = min(pre),
            max_PRE = max(pre),
            med_sd_APE_within = median(sd_ape)) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  ungroup()

error_metrics$sc3 <- as.factor(error_metrics$sc3)
error_metrics$sc3 <- factor(error_metrics$sc3, levels(error_metrics$sc3)[c(6, 5,4,1,3,2,8,7,11,12,9,10)])

summary(error_metrics$N)
summary(error_metrics$nrmse)

# Figure S10: accuracy metrics plots for the females
f_ac_plot <- ggplot(error_metrics %>%
                      filter(data == "Females") %>%
                      dplyr::select(sc3, et, med_PRE, nrmse) %>%
                      rename(`Median PRE` = med_PRE,
                             NRMSE = nrmse) %>%
                      gather(metric, value, -et, -sc3)) +
  geom_col(aes(x = sc3, y = value, fill = et), color = "black",
           position = position_dodge2()) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  guides(fill = guide_legend(title = "Parameter Estimation")) + 
  xlab("Scenario") + ylab("Metric Value")

# Figure S11: accuracy metrics plots for the males
m_ac_plot <- ggplot(error_metrics %>%
                      filter(data == "Males") %>%
                      dplyr::select(sc3, et, med_PRE, nrmse) %>%
                      rename(`Median PRE` = med_PRE,
                             NRMSE = nrmse) %>%
                      gather(metric, value, -et, -sc3)) +
  geom_col(aes(x = sc3, y = value, fill = et), color = "black",
           position = position_dodge2()) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  guides(fill = guide_legend(title = "Parameter Estimation")) + 
  xlab("Scenario") + ylab("Metric Value")



# ------------------------------------
# error metrics table

sr_error_metrics <- error_metrics_scens %>%
  dplyr::select(eofsr, ped_rep, mating, nsampyrs, ssf, ckmr_ssmult,
                et, year, hze, hzt, nll_type, fec_known, ckmr_seed, 
                data, sim_value, est_value, id, scenario, sc3) %>%
  mutate(data2 = data) %>%
  spread(data, est_value) %>%
  rename(fem_est = Females, male_est = Males) %>%
  spread(data2, sim_value) %>%
  rename(fem_sim = Females, male_sim = Males) %>%
  group_by(et, year, id, sc3) %>%
  summarize(fem_est = sum(fem_est, na.rm = TRUE),
            male_est = sum(male_est, na.rm = TRUE),
            fem_sim = sum(fem_sim, na.rm = TRUE),
            male_sim = sum(male_sim, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sim_sr = male_sim/(male_sim+fem_sim)*100,
         est_sr = male_est/(male_est+fem_est)*100,
         pre = get_perror(sim = sim_sr, est = est_sr),
         ape = abs(est_sr - sim_sr)/sim_sr * 100) %>%
  group_by(sc3, et, id) %>%
  summarize(pre = median(pre),
            ape = median(ape),
            sd_ape = sd(ape),
            med_sim = median(sim_sr),
            med_est = median(est_sr),
            squared_diff = (med_est-med_sim)^2) %>%
  ungroup() %>%
  group_by(sc3, et) %>%
  summarize(med_PRE = median(pre),
            med_APE = median(ape),
            iqr_PRE = IQR(pre),
            N = length(unique(id)),
            meanobs = mean(med_sim),
            rmse = sqrt((1/N) * sum(squared_diff)),
            nrmse = rmse/meanobs,
            min_PRE = min(pre),
            max_PRE = max(pre),
            med_sd_APE_within = median(sd_ape)) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  ungroup() %>%
  mutate(data = "SexRatio",
         .after = et)

sr_error_metrics$sc3 <- as.factor(sr_error_metrics$sc3)
sr_error_metrics$sc3 <- factor(sr_error_metrics$sc3, levels(sr_error_metrics$sc3)[c(6, 5,4,1,3,2,8,7,11,12,9,10)])

names(sr_error_metrics) == names(error_metrics)

error_metrics_full <- error_metrics %>%
  bind_rows(sr_error_metrics) %>%
  pivot_wider(id_cols = c(sc3, et), 
              names_from = data, 
              values_from = c("N", "med_PRE", "iqr_PRE", "nrmse")) %>%
  dplyr::select(-N_Males, -N_SexRatio) %>%
  mutate(med_PRE_Males = round(med_PRE_Males,2),
         med_PRE_Females = round(med_PRE_Females,2),
         med_PRE_SexRatio = round(med_PRE_SexRatio,2),
         iqr_PRE_Females = round(iqr_PRE_Females,1),
         iqr_PRE_Males = round(iqr_PRE_Males,1),
         iqr_PRE_SexRatio = round(iqr_PRE_SexRatio,1),) %>%
  rename(N = N_Females,
         Scenario = sc3,
         `Parameter Estimation` = et) %>%
  relocate(Scenario, `Parameter Estimation`, N, med_PRE_Females, iqr_PRE_Females, 
           nrmse_Females, med_PRE_Males, iqr_PRE_Males, nrmse_Males,
           med_PRE_SexRatio, iqr_PRE_SexRatio, nrmse_SexRatio) %>%
  arrange(Scenario, `Parameter Estimation`)

names(error_metrics_full)[4:12] <- rep(c("PRE", "IQR", "NRMSE"), 3)



#--------------------------------
# Figures A12-A14: PRE by simulated pedigree

# female data

fdat <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult,
         !(om_sc == "skipped spawning" & em_sc == "HSPs only"),
         data == "Females",
  ) %>%
  mutate(sc3 = case_when(
    em_sc != "Base" ~ em_sc,
    om_sc == "Base" ~ "Base",
    TRUE ~ om_sc
  )) %>%
  group_by(data, ped_rep, ckmr_seed, et, sc3, id) %>%
  summarize(perror = median(perror),
            Nyears = n(),
            rmse = sqrt((1/Nyears) * sum(squared_diff))) %>%
  ungroup() %>%
  mutate(ped_rep = as.factor(ped_rep),
         #perror = scale(perror),
         #rmse = scale(rmse)
  ) %>%
  filter(sc3 != "HSPs only") # nothing interesting here, so don't need to show this
dim(fdat)
summary(fdat$Nyears)

fdat$sc3 <- factor(fdat$sc3,
                   levels = c("Base",
                              "skipped spawning",
                              "mRV flatter",
                              "mRV steeper",
                              "est Hmt < sim Hmt",
                              "sim Hmt < est Hmt"))

#table(fdat$ped_rep, fdat$et, fdat$sc3)

fdat <- fdat %>% group_by(ped_rep, et, sc3) %>% 
  mutate(N = n())

fdat <- fdat %>% filter(N > 3)
fdat <- fdat %>% droplevels()


# male data

mdat <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult,
         !(om_sc == "skipped spawning" & em_sc == "HSPs only"),
         data == "Males",
  ) %>%
  mutate(sc3 = case_when(
    em_sc != "Base" ~ em_sc,
    om_sc == "Base" ~ "Base",
    TRUE ~ om_sc
  )) %>%
  group_by(data, ped_rep, ckmr_seed, et, sc3, id) %>%
  summarize(perror = median(perror),
            Nyears = n(),
            rmse = sqrt((1/Nyears) * sum(squared_diff))) %>%
  ungroup() %>%
  mutate(ped_rep = as.factor(ped_rep),
         #perror = scale(perror),
         #rmse = scale(rmse)
  ) %>%
  filter(sc3 != "HSPs only") # nothing interesting here, so don't need to show this

mdat$sc3 <- factor(mdat$sc3,
                   levels = c("Base",
                              "skipped spawning",
                              "mRV flatter",
                              "mRV steeper",
                              "est Hmt < sim Hmt",
                              "sim Hmt < est Hmt"))

mdat <- mdat %>% group_by(ped_rep, et, sc3) %>% 
  mutate(N = n())

mdat <- mdat %>% filter(N > 3)
mdat <- mdat %>% droplevels()

### Figure S12: female plot

fpr_plot <- ggplot(fdat %>% filter(perror < 1500), # omits two points for est Hmt < sim Hmt
                   aes(x = ped_rep, y = perror, fill = et)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.2,
               size = 0.1,
               position = position_dodge(preserve = "single")) +
  coord_flip() +
  labs(#title = "Female Abundance",
    x = "Simulated Population", y = "Percent Relative Error") +
  theme_minimal() +
  facet_wrap(~sc3, scales = "free_x", ncol = 6,
             axis.labels = "margins") +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(4, "mm"),
        axis.text = element_text(size = 6)) +
  guides(fill = guide_legend(title = "EM Scenario"))


### Figure S13: male plot

mpr_plot <- ggplot(mdat %>% filter(perror < 1500), # omits two points for est Hmt < sim Hmt
                   aes(x = ped_rep, y = perror, fill = et)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.2,
               size = 0.1,
               position = position_dodge(preserve = "single")) +
  coord_flip() +
  labs(#title = "Male Abundance",
    x = "Simulated Population", y = "Percent Relative Error") +
  theme_minimal() +
  facet_wrap(~sc3, scales = "free_x", ncol = 6,
             axis.labels = "margins") +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(4, "mm"),
        axis.text = element_text(size = 6)) +
  guides(fill = guide_legend(title = "EM Scenario"))


### Figure S14: Marginal means plot for scenario effect

library(glmmTMB)
library(MuMIn)
library(emmeans)

### Females


# fit1 <- glmmTMB(perror ~ et,
#                 data = fdat, family = gaussian())
# fit2 <- glmmTMB(perror ~ et + sc3,
#                 data = fdat, family = gaussian())
# fit3 <- glmmTMB(perror ~ et + sc3 + (1|ped_rep),
#                 data = fdat, family = gaussian())
# fit4 <- glmmTMB(perror ~ et * sc3,
#                 data = fdat, family = gaussian())
fit5 <- glmmTMB(perror ~ et * sc3 + (1|ped_rep),
                data = fdat, family = gaussian())

# model.sel(fit1, fit2, fit3, fit4, fit5) # fit5 is best

# VarCorr(fit5)
# summary(fit5)
# residuals <- residuals(fit5)
# plot(residuals)
# confint(fit5)


### visualize fixed effects
fixef_df <- as.data.frame(summary(fit5)$coefficients$cond) %>%
  mutate(CI_lower = Estimate - 1.96 * `Std. Error`,
         CI_upper = Estimate + 1.96 * `Std. Error`,
         term = rownames(.))

# ggplot(fixef_df, aes(x = reorder(term, Estimate), y = Estimate)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper)) +
#   coord_flip() + 
#   labs(title = "Fixed Effects Coefficients for perror",
#        x = "Effect", y = "Estimate") +
#   theme_minimal()


### visualize random effects
ranef_vals <- ranef(fit5)[[1]]$ped_rep

# Convert to data frame for plotting
ranef_df <- data.frame(
  ped_rep = sort(unique(fdat$ped_rep)),
  random_intercept = ranef_vals[, 1]
)

# Reorder ped_rep based on the values of random_intercept
ranef_df$ped_rep <- factor(ranef_df$ped_rep, levels = ranef_df$ped_rep[order(ranef_df$random_intercept)])

# Plot random effects
# ggplot(ranef_df, aes(x = ped_rep, y = random_intercept)) +
#   geom_point() +
#   coord_flip() +
#   labs(title = "Random Effects: Intercepts with Confidence Intervals for Each ped_rep",
#        x = "ped_rep (Population)", y = "Random Intercept") +
#   theme_minimal()

# marginal means for females
femm <- plot(emmeans(fit5, ~ sc3 * et), comparisons = TRUE, plotit = F)


### Males

# fit6 <- glmmTMB(perror ~ et,
#                 data = mdat, family = gaussian())
# fit7 <- glmmTMB(perror ~ et + sc3,
#                 data = mdat, family = gaussian())
# fit8 <- glmmTMB(perror ~ et + sc3 + (1|ped_rep),
#                 data = mdat, family = gaussian())
fit9 <- glmmTMB(perror ~ et * sc3,
                data = mdat, family = gaussian())
# fit10 <- glmmTMB(perror ~ et * sc3 + (1|ped_rep),
#                 data = mdat, family = gaussian())
# 
# model.sel(fit6, fit7, fit8, fit9, fit10) # fit9 is best, no random effects!

# marginal means for males
memm <- plot(emmeans(fit9, ~ sc3 * et), comparisons = TRUE, plotit = F)

full_emm <- femm %>% mutate(sex = "Female") %>%
  bind_rows(memm %>% mutate(sex = "Male"))

emmeans_plot <- ggplot(full_emm, aes(x = et, y = the.emmean, color = sc3, group = sc3)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  theme_minimal() +  # Use a minimal theme
  facet_wrap(~sex) +
  labs(
    x = "Parameter Estimation",
    y = "Estimated Marginal Mean PRE",
    color = "Scenario") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))




#--------------------------------
# Figure S15: number sampled and kin pairs found when
# varying the proportion of males in the CKMR sample

samp_sum_suppl_plot <- ggplot(samp_summary %>%
                                filter(LifeStage != "Juv") %>%
                                mutate(Sex = ifelse(Sex == "F", "Female", "Male")) %>%
                                filter(om_sc %in% c("Base","skipped spawning"), em_sc == "Base", et == base_et,
                                       nsampyrs == base_ckmr_nsampyrs, ckmr_ssmult == base_ckmr_ssmult)) +
  geom_boxplot(aes(x = om_sc, y = N, fill = sample_by_sex),
               position = position_dodge(preserve = "single"),
               outlier.size = 1.4, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.2,
               size = 0.4) +
  facet_wrap(~Sex, scales = "free_y") +
  xlab("") + ylab("Number Sampled") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                    labels = c("Base", "More males", "No males")) +
  theme(legend.position = "none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  labs(tag = "A")

ckmr_pairs_suppl <- ckmr_obs_samps_long %>% 
  filter(em_sc == "Base", nsampyrs == base_ckmr_nsampyrs, om_sc %in% c("Base","skipped spawning"), 
         ckmr_ssmult == base_ckmr_ssmult, et == base_et)

ckmr_pairs_suppl_TEX <- ckmr_pairs_suppl
levels(ckmr_pairs_suppl_TEX$pair_type) <- c(
  POD = TeX("POP$^\\neq$"),
  POS = TeX("POP$^=$"),
  HSD = TeX("HSP$^\\neq$"),
  HSS = TeX("HSP$^=$"),
  Unrelated = TeX("Unrelated")
)

pairs_plot_suppl_plot <- ggplot(ckmr_pairs_suppl_TEX) +
  geom_boxplot(aes(x = om_sc, y = count, fill = sample_by_sex),
               position = position_dodge(preserve = "single"),
               outlier.size = 1.4, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.2,
               size = 0.4) +
  facet_wrap(~pair_type, labeller = label_parsed,
             scales = "free_y", ncol = 3) +
  ylab("Number of Pairs") + xlab("") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                    labels = c("Base", "More males", "No males")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        strip.text.x = element_text(size = 12)) +
  guides(fill = guide_legend(title = "Male Sampling")) +
  labs(tag = "B")


male_samps_suppl_plot1 <- samp_sum_suppl_plot + pairs_plot_suppl_plot + 
  plot_layout(
    nrow = 2,
    heights = c(0.5,1)) &
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text = element_text(size = 10))



#--------------------------------
# Figure S16: PRE for varying the proportion of males in the CKMR sample

samps_bysex <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, ckmr_ssmult == base_ckmr_ssmult,
         om_sc %in% c("Base","skipped spawning"),
         em_sc == "Base") %>%
  group_by(om_sc, sample_by_sex, et, data, id) %>%
  summarize(perror = median(perror))


sbys_fem_p <- ggplot(samps_bysex %>% filter(data == "Females"), aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = sample_by_sex),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               size = 0.4) +
  facet_wrap(~om_sc) + 
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                    labels = c("Base", "More males", "No males")) +
  xlab("") + ylab("Percent Relative Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  ggtitle("Females") +
  labs(tag = "A")


sbys_mal_p <- ggplot(samps_bysex %>% filter(data == "Males"), aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = sample_by_sex),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               size = 0.4) +
  facet_wrap(~om_sc) + 
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF"),
                    labels = c("Base", "More males", "No males")) +
  xlab("") + ylab("Percent Relative Error") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  ggtitle("Males") +
  labs(tag = "B")

male_samps_suppl_plot2 <- sbys_fem_p + sbys_mal_p + 
  plot_layout(axis_titles = "collect",
              guides = "collect",
              nrow = 1) &
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text = element_text(size = 10)) &
  guides(fill = guide_legend(title = "Male Sampling"))

