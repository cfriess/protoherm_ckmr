
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
em_sc_levels <- c("Base", "est Hmt < sim Hmt", "sim Hmt < est Hmt")

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
sim_pars$scenario <- factor(sim_pars$scenario, levels(sim_pars$scenario)[c(1,5,2,6,3,4,7,8)])

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
                                size = 12, 
                                fontface = "plain", 
                                angle = 90)

xtitle <- ggdraw() + draw_label("CKMR Sample Size", 
                                size = 12, 
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
# Figure S3: PRE by iteration for males

# note filtering for only those where all years are < 1000 PRE
# results in omitting 10 iterations (all for eHmt < sHmt). 
# no interesting patterns, just really high errors:

# abund_res %>% 
#   filter(scen_type == "Main", data == "Females") %>%
#   group_by(scenario, et, id) %>%
#   mutate(tmp_id = cur_group_id()) %>%  # first get a unique group number per scenario+id combo
#   group_by(scenario) %>%
#   mutate(tmp_id = dense_rank(tmp_id)) %>% # then re-rank within scenario so it starts at 1
#   ungroup() %>%
#   group_by(scenario, et, id) %>%
#   filter(any(perror > 1000)) %>%
#   ungroup() %>%
#   distinct(scenario, et, id)

pre_by_year_females <- abund_res %>% 
  filter(scen_type == "Main", data == "Females") %>%
  group_by(scenario, et, id) %>%
  mutate(tmp_id = cur_group_id()) %>%  # first get a unique group number per scenario+id combo
  group_by(scenario) %>%
  mutate(tmp_id = dense_rank(tmp_id)) %>% # then re-rank within scenario so it starts at 1
  ungroup() %>%
  group_by(scenario, et, id) %>%
  filter(all(perror < 1000)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = perror, color = as.factor(tmp_id))) +
  geom_line() +
  facet_grid(et~scenario) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "PRE", x = "Year") +
  theme(legend.position =  "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 9))


#--------------------------------
# Figure S4: PRE by iteration for males

# note filtering for only those where all years are < 500 PRE
# results in omitting 20 iterations (6 for eHmt < sHmt fixed trans & est mRV
# and 16 for mRV steeper est trans & mRV). no interesting patterns, just
# really high errors:

# abund_res %>% 
#   filter(scen_type == "Main", data == "Males") %>%
#   group_by(scenario, et, id) %>%
#   mutate(tmp_id = cur_group_id()) %>% 
#   group_by(scenario) %>%
#   mutate(tmp_id = dense_rank(tmp_id)) %>% 
#   ungroup() %>%
#   group_by(scenario, et, id) %>%
#   filter(any(perror > 500)) %>%
#   ungroup() %>%
#   distinct(scenario, et, id)

pre_by_year_males <- abund_res %>% 
  filter(scen_type == "Main", data == "Males") %>%
  group_by(scenario, et, id) %>%
  mutate(tmp_id = cur_group_id()) %>%  # first get a unique group number per scenario+id combo
  group_by(scenario) %>%
  mutate(tmp_id = dense_rank(tmp_id)) %>% # then re-rank within scenario so it starts at 1
  ungroup() %>%
  group_by(scenario, et, id) %>%
  filter(all(perror < 500)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = perror, color = as.factor(tmp_id))) +
  geom_line() +
  facet_grid(et~scenario) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "PRE", x = "Year") +
  theme(legend.position =  "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 9))


#--------------------------------
# Figure S5: distribution of PREs

dat <- abund_res %>% filter(scen_type == "Main") %>%
  group_by(id, ped_rep, scenario, et, data) %>%
  summarize(perror = median(perror), .groups = "drop")

plot_settings <- list(
  geom_vline(xintercept = 0, linetype = "dashed"),
  facet_wrap(~et, scales = "free"),
  theme(plot.title = element_text(size = 9),
        strip.text = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.title = element_blank())  # remove individual axis labels
)

p1 <- dat %>% filter(data == "Males", scenario == "Base") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#F8766D") + plot_settings + ggtitle("Base")

p2 <- dat %>% filter(data == "Males", scenario == "sHmt < eHmt") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#7CAE00") + plot_settings + ggtitle("sHmt < eHmt")

p3 <- dat %>% filter(data == "Males", scenario == "eHmt < sHmt") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#00BFC4") + plot_settings + ggtitle("eHmt < sHmt")

p4 <- dat %>% filter(data == "Males", scenario == "skipped spawning") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#C77CFF") + plot_settings + ggtitle("skipped spawning")

p5 <- dat %>% filter(data == "Males", scenario == "mRV flatter") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#E58700") + plot_settings + ggtitle("mRV flatter")

p6 <- dat %>% filter(data == "Males", scenario == "mRV steeper") %>%
  ggplot(aes(x = perror)) + geom_histogram(fill = "#00BE67") + plot_settings + ggtitle("mRV steeper")


#--------------------------------
# Figure S6: Number of individuals sampled for CKMR

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
# Figure S7: Number of observed kin pairs by CKMR sampling

ckmr_pairs_by_sampling <- ckmr_samps_long %>% 
  filter(scenario %in% c("Base", "Suppl Sampling Intensity"))

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
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size = 12)
  ) +
  xlab("CKMR Sample Size")


#--------------------------------
# Figure S8: ckmr pairs by EM scenario


ckmr_pairs_by_scenario2 <- ckmr_obs_samps_long %>% 
  filter(scenario %in% c("Base", "sHmt < eHmt", "eHmt < sHmt"))

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
  guides(fill = guide_legend(title = "Scenario"))



#--------------------------------
# Figure S9: estimated parameters by EM scenario

layout3 <- "
BCF
DEA
" 

pars_nll <- sim_pars %>%
  filter(scenario %in% c("Base", "sHmt < eHmt", "eHmt < sHmt")) 

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
                   outlier.shape = 21, outlier.stroke = 0.3,
                   linewidth = 0.2) +
      scale_fill_discrete(
        labels = c(
          "Base" = "Base",
          "est Hmt < sim Hmt" = expression(H[mt]^{"EM"} < H[mt]^{"OM"}),
          "sim Hmt < est Hmt" = expression(H[mt]^{"OM"} < H[mt]^{"EM"})
        )
      ) +
      facet_nested(. ~ par_name_long) +
      scale_x_discrete(labels = scales::wrap_format(15)) +
      ylab("Percent Relative Error") +
      guides(fill = guide_legend(title = "Scenario")) 
    
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
# Figure S10: estimated parameters by OM scenario

pars_mRV <- sim_pars %>%
  filter(scenario %in% c("Base", "skipped spawning", "mRV flatter", "mRV steeper")
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
                   outlier.shape = 21, outlier.stroke = 0.3,
                   linewidth = 0.2) +
      scale_fill_discrete(
        labels = om_labels
      ) +
      facet_nested(. ~ par_name_long) +
      scale_x_discrete(labels = scales::wrap_format(15)) +
      ylab("Percent Relative Error") +
      guides(fill = guide_legend(title = "Scenario")) 
    
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
# Figure S11: parameter correlations

cor_dat <- cor_res_et %>%
  filter(scenario %in% c("Base", "Suppl Sampling Intensity")) %>%
  dplyr::select(nsampyrs, ckmr_ssmult, rbar_tsd_cor:tsd_exp_cor,fec_known, om_sc) %>%
  gather(comp, cor, rbar_tsd_cor:tsd_exp_cor) %>%
  mutate(comp = case_when(
    comp == "rbar_tsd_cor" ~ "Recruitment & \nTransition Fct",
    comp == "rbar_exp_cor" ~ "Recruitment & \nMale RV",
    comp == "tsd_exp_cor" ~ "Transition Fct & \nMale RV"
  ))


cor_plot <- ggplot(cor_dat) +
  geom_boxplot(aes(x = ckmr_ssmult, y = cor), fill = "grey") +
  facet_wrap(~comp) +
  scale_fill_discrete(
    labels = om_labels
  ) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(size = 11)) +
  xlab("CKMR Sample Size") +
  ylab("Correlation")



#--------------------------------
# Figures A12-A13: PRE and NRMSE comparison

error_metrics_scens <- abund_res %>%
              filter(scen_type == "Main" | scenario == "Suppl Sampling Intensity"#,
                     # !(ckmr_ssmult == base_ckmr_ssmult &
                     #     nsampyrs== base_ckmr_nsampyrs)
                     ) %>%
              mutate(scenario = case_when(
                nsampyrs == "03 Yrs" & ckmr_ssmult == "50%" ~ "3 Yrs, 50% BaseN",
                nsampyrs == "03 Yrs" & ckmr_ssmult == "150%" ~ "3 Yrs, 150% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "50%" ~ "10 Yrs, 50% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "100%" ~ "10 Yrs, 100% BaseN",
                nsampyrs == "10 Yrs" & ckmr_ssmult == "150%" ~ "10 Yrs, 150% BaseN",
                TRUE ~ scenario
              ))

error_metrics <- error_metrics_scens %>%
  group_by(scenario, et, data, id, ped_rep, ckmr_seed) %>% # summarize by id
  summarize(pre = median(perror),
            ape = median(aperror),
            sd_ape = sd(aperror),
            med_sim = median(sim_value),
            med_est = median(est_value),
            squared_diff = (med_est-med_sim)^2) %>%
  ungroup() %>%
  group_by(scenario, et, data) %>% # summarize by scenario
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

error_metrics$scenario <- as.factor(error_metrics$scenario)
error_metrics$scenario <- factor(error_metrics$scenario, levels(error_metrics$scenario)[c(6,10,7,11,8,9,5,4,3,1,2)])

summary(error_metrics$N)
summary(error_metrics$nrmse)



# Figure S12: accuracy metrics plots for the females
f_ac_plot <- ggplot(error_metrics %>%
                      filter(data == "Females") %>%
                      dplyr::select(scenario, et, med_PRE, nrmse) %>%
                      rename(`Median PRE` = med_PRE,
                             NRMSE = nrmse) %>%
                      gather(metric, value, -et, -scenario)) +
  geom_col(aes(x = scenario, y = value, fill = et), color = "black",
           position = position_dodge2()) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(title = "Parameter Estimation", nrow = 2, byrow = T)) + 
  xlab("Scenario") + ylab("Metric Value")

# Figure S13: accuracy metrics plots for the males
m_ac_plot <- ggplot(error_metrics %>%
                      filter(data == "Males") %>%
                      dplyr::select(scenario, et, med_PRE, nrmse) %>%
                      rename(`Median PRE` = med_PRE,
                             NRMSE = nrmse) %>%
                      gather(metric, value, -et, -scenario)) +
  geom_col(aes(x = scenario, y = value, fill = et), color = "black",
           position = position_dodge2()) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = "bottom",
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(title = "Parameter Estimation", nrow = 2, byrow = T)) + 
  xlab("Scenario") + ylab("Metric Value")




#--------------------------------
### Figure S14: ICC plot

library(lme4)

dat <- abund_res %>% filter(scen_type == "Main") %>%
  group_by(id, ped_rep, scenario, et, data) %>%
  summarize(perror = median(perror), .groups = "drop")

icc_results <- dat %>%
  group_by(scenario, et, data) %>%
  group_modify(~ {
    m <- lmer(perror ~ 1 + (1|ped_rep), data = .x)
    vc <- as.data.frame(VarCorr(m))
    icc <- vc$vcov[1] / sum(vc$vcov)
    tibble(icc = icc, singular = isSingular(m))
  })

icc_results$scenario <- as.factor(icc_results$scenario)
icc_results$scenario <- factor(icc_results$scenario, levels = 
                                 levels(icc_results$scenario)[c(1,5,2,6,3,4)])

icc_plot <- ggplot(icc_results) +
  geom_col(aes(x = scenario, y = icc, fill = et),
           position = position_dodge2()) +
  facet_wrap(~data) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  guides(fill = guide_legend(title = "Parameter Estimation",
                             nrow = 2, byrow = T)) +
  xlab("Scenario") + ylab("Intraclass Correlation Coefficient")


#--------------------------------
# Figures S15-S16: PRE by simulated pedigree

# female data

fdat <- abund_res %>%
  filter(scen_type == "Main",
         data == "Females",
  ) %>%
  group_by(data, ped_rep, ckmr_seed, et, scenario, id) %>%
  summarize(perror = median(perror),
            Nyears = n(),
            rmse = sqrt((1/Nyears) * sum(squared_diff))) %>%
  ungroup() %>%
  mutate(ped_rep = as.factor(ped_rep),
         #perror = scale(perror),
         #rmse = scale(rmse)
  ) 
dim(fdat)
summary(fdat$Nyears)

fdat$scenario <- factor(fdat$scenario,
                        levels = c("Base",
                                   "sHmt < eHmt",
                                   "eHmt < sHmt",
                                   "skipped spawning",
                                   "mRV flatter",
                                   "mRV steeper"
                        ))


fdat <- fdat %>% group_by(ped_rep, et, scenario) %>% 
  mutate(N = n())

fdat <- fdat %>% filter(N > 3)
fdat <- fdat %>% droplevels()


# male data

mdat <- abund_res %>%
  filter(scen_type == "Main",
         data == "Males",
  ) %>%
  group_by(data, ped_rep, ckmr_seed, et, scenario, id) %>%
  summarize(perror = median(perror),
            Nyears = n(),
            rmse = sqrt((1/Nyears) * sum(squared_diff))) %>%
  ungroup() %>%
  mutate(ped_rep = as.factor(ped_rep),
         #perror = scale(perror),
         #rmse = scale(rmse)
  )

mdat$scenario <- factor(mdat$scenario,
                        levels = c("Base",
                                   "sHmt < eHmt",
                                   "eHmt < sHmt",
                                   "skipped spawning",
                                   "mRV flatter",
                                   "mRV steeper"
                        ))

mdat <- mdat %>% group_by(ped_rep, et, scenario) %>% 
  mutate(N = n())

mdat <- mdat %>% filter(N > 3)
mdat <- mdat %>% droplevels()

### Figure S15: female plot

fpr_plot <- ggplot(fdat %>% filter(perror < 1500,
                                   ped_rep %in% 1:20), # omits one point for est Hmt < sim Hmt
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
  facet_wrap(~scenario, scales = "free_x", ncol = 6,
             axis.labels = "margins") +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(4, "mm"),
        axis.text = element_text(size = 6)) +
  guides(fill = guide_legend(title = "EM Scenario"))


### Figure S16: male plot

mpr_plot <- ggplot(mdat %>% filter(ped_rep %in% 1:20),
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
  facet_wrap(~scenario, scales = "free_x", ncol = 6,
             axis.labels = "margins") +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(4, "mm"),
        axis.text = element_text(size = 6)) +
  guides(fill = guide_legend(title = "EM Scenario"))


#--------------------------------
# Figure S17: PRE by number of CKMR samples and number of sampling yrs


err_ckmr_sampling <- abund_res %>%
  filter(scenario %in% c("Base","Suppl Sampling Intensity")) %>%
  group_by(data, et, ckmr_ssmult, nsampyrs, id) %>%
  summarize(perror = median(perror)) %>%
  ungroup()


sp1 <- err_ckmr_sampling %>% filter(et == "fixed trans & mRV") %>%
  ggplot(aes(x = ckmr_ssmult, y = perror, fill = nsampyrs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.4,
               linewidth = 0.3) +
  facet_wrap( ~ data, scales = "fixed") +
  xlab("CKMR Sample Size") +
  ylab("Percent Relative Error") +
  ggtitle(expression(bold("fixed" ~ bolditalic(P)^{"F" %->% "M"} ~ " & RV"["M"]))) +
  #ggtitle("fixed trans & mRV") +
  coord_cartesian(ylim = c(-55,150)) +
  labs(tag = "A") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

sp2 <- err_ckmr_sampling %>% filter(et == "est trans & fixed mRV") %>%
  ggplot(aes(x = ckmr_ssmult, y = perror, fill = nsampyrs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.4,
               linewidth = 0.3) +
  facet_wrap( ~ data, scales = "fixed") +
  xlab("CKMR Sample Size") +
  ylab("Percent Relative Error") +
  ggtitle(expression(bold("free" ~ bolditalic(P)^{"F" %->% "M"} ~ ", fixed RV"["M"]))) +
  #ggtitle("est trans & fixed mRV") +
  coord_cartesian(ylim = c(-55,150)) + # omits 936, highest outliers were for male 10 yrs 50% sampling est both; max 1129.9773
  labs(tag = "B") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

sp3 <- err_ckmr_sampling %>% filter(et == "fixed trans & est mRV") %>%
  ggplot(aes(x = ckmr_ssmult, y = perror, fill = nsampyrs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.4,
               linewidth = 0.3) +
  facet_wrap( ~ data, scales = "fixed") +
  xlab("CKMR Sample Size") +
  ylab("Percent Relative Error") +
  ggtitle(expression(bold("fixed" ~ bolditalic(P)^{"F" %->% "M"} ~ ", free RV"["M"]))) +
  #ggtitle("fixed trans & est mRV") +
  coord_cartesian(ylim = c(-55,150)) + # omits 936, highest outliers were for male 10 yrs 50% sampling est both; max 1129.9773
  labs(tag = "C") +
  theme(plot.title = element_text(hjust = 0.5))

sp4 <- err_ckmr_sampling %>% filter(et == "est trans & mRV") %>%
  ggplot(aes(x = ckmr_ssmult, y = perror, fill = nsampyrs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.4,
               linewidth = 0.3) +
  facet_wrap( ~ data, scales = "fixed") +
  xlab("CKMR Sample Size") +
  ylab("Percent Relative Error") +
  ggtitle(expression(bold("free" ~ bolditalic(P)^{"F" %->% "M"} ~ " & RV"["M"]))) +
  #ggtitle("est trans & mRV") +
  coord_cartesian(ylim = c(-55,150)) + # omits 936, highest outliers were for male 10 yrs 50% sampling est both; max 1129.9773
  labs(tag = "D") +
  theme(plot.title = element_text(hjust = 0.5))

layout2 <- "
ABC
ADE
FFF
"

ckmr_sampling_plot <- ytitle + sp1 + sp2 + sp3 + sp4 + xtitle +
  plot_layout(
    guides = "collect",
    design = layout2,
    widths = c(0.05,0.475,0.475),
    heights = c(0.475,0.475,0.05)) &
  theme(legend.position = "none",
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.title = element_blank(),
        axis.text = element_text(size = 10))



#--------------------------------
# Figure S18: number sampled and kin pairs found when
# varying the proportion of males in the CKMR sample

samp_sum_suppl_plot <- ggplot(samp_summary %>%
                                filter(LifeStage != "Juv") %>%
                                mutate(Sex = ifelse(Sex == "F", "Female", "Male")) %>%
                                filter(scenario %in% c("Base","skipped spawning","Suppl Scen Male & Female Sampling"))) +
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
  filter(scenario %in% c("Base","skipped spawning","Suppl Scen Male & Female Sampling"))

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
# Figure S19: PRE for varying the proportion of males in the CKMR sample

samps_bysex <- abund_res %>%
  filter(scenario %in% c("Base","skipped spawning","Suppl Scen Male & Female Sampling")) %>%
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




# ------------------------------------
# error metrics table

sr_error_metrics <- error_metrics_scens %>%
  dplyr::select(eofsr, ped_rep, mating, nsampyrs, ssf, ckmr_ssmult,
                et, year, hze, hzt, nll_type, fec_known, ckmr_seed, 
                data, sim_value, est_value, id, scenario, scenario) %>%
  mutate(data2 = data) %>%
  spread(data, est_value) %>%
  rename(fem_est = Females, male_est = Males) %>%
  spread(data2, sim_value) %>%
  rename(fem_sim = Females, male_sim = Males) %>%
  group_by(et, year, id, scenario) %>%
  summarize(fem_est = sum(fem_est, na.rm = TRUE),
            male_est = sum(male_est, na.rm = TRUE),
            fem_sim = sum(fem_sim, na.rm = TRUE),
            male_sim = sum(male_sim, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sim_sr = male_sim/(male_sim+fem_sim)*100,
         est_sr = male_est/(male_est+fem_est)*100,
         pre = get_perror(sim = sim_sr, est = est_sr),
         ape = abs(est_sr - sim_sr)/sim_sr * 100) %>%
  group_by(scenario, et, id) %>%
  summarize(pre = median(pre),
            ape = median(ape),
            sd_ape = sd(ape),
            med_sim = median(sim_sr),
            med_est = median(est_sr),
            squared_diff = (med_est-med_sim)^2) %>%
  ungroup() %>%
  group_by(scenario, et) %>%
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

sr_error_metrics$scenario <- as.factor(sr_error_metrics$scenario)
sr_error_metrics$scenario <- factor(sr_error_metrics$scenario, 
                                    levels(sr_error_metrics$scenario)[c(6,10,7,11,8,9,5,4,3,1,2)])

names(sr_error_metrics) == names(error_metrics)

error_metrics_full <- error_metrics %>%
  bind_rows(sr_error_metrics) %>%
  pivot_wider(id_cols = c(scenario, et), 
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
         Scenario = scenario,
         `Parameter Estimation` = et) %>%
  relocate(Scenario, `Parameter Estimation`, N, med_PRE_Females, iqr_PRE_Females, 
           nrmse_Females, med_PRE_Males, iqr_PRE_Males, nrmse_Males,
           med_PRE_SexRatio, iqr_PRE_SexRatio, nrmse_SexRatio) %>%
  arrange(Scenario, `Parameter Estimation`)

names(error_metrics_full)[4:12] <- rep(c("PRE", "IQR", "NRMSE"), 3)

#--------------------------------
# CKMR sampling table

# n for each of these is 400
obs_summary <- ckmr_obs_samps_long %>% 
  filter(sample_by_sex == samp_sex_base, et == base_et) %>%
  group_by(scenario, nsampyrs, ckmr_ssmult, pair_type) %>%
  summarize(median_count = round(median(count),0),
            #min_count = min(count),
            #max_count = max(count),
            iqr_count = round(IQR(count),0)#,
            #n = n()
  ) %>%
  mutate(count = paste0(format(median_count, big.mark = ",", scientific = FALSE), " (", 
                        format(iqr_count, big.mark = ",", scientific = FALSE, trim = TRUE), ")")) %>%
  dplyr::select(-median_count, -iqr_count) %>%
  pivot_wider(id_cols = c(scenario, nsampyrs, ckmr_ssmult), 
              names_from = pair_type, 
              values_from = c("count"))

# n for each of these is 400
sample_summary <- samp_summary %>%
  filter(sample_by_sex == samp_sex_base, et == base_et, LifeStage == "Adult") %>%
  group_by(scenario, nsampyrs, ckmr_ssmult, Sex) %>%
  summarize(median_N = round(median(N),0),
            #min_N = min(N),
            #max_N = max(N),
            iqr_N = round(IQR(N),0)#,
            #n = n()
  ) %>%
  mutate(count = paste0(format(median_N, big.mark = ",", scientific = FALSE), " (", 
                        format(iqr_N, big.mark = ",", scientific = FALSE, trim = TRUE), ")")) %>%
  dplyr::select(-median_N, -iqr_N) %>%
  pivot_wider(id_cols = c(scenario, nsampyrs, ckmr_ssmult), 
              names_from = Sex, 
              values_from = c("count"))


samps_table_data <- sample_summary %>%
  left_join(obs_summary, 
            by = c("scenario", "nsampyrs", "ckmr_ssmult")
  ) %>%
  relocate(scenario, nsampyrs, ckmr_ssmult, 
           F, M, POS, POD, HSS, HSD, Unrelated
  ) %>%
  mutate(scenario = ifelse(scenario == "Suppl Sampling Intensity",
         "Varied Sampling", scenario))

samps_table_data$scenario <- as.factor(samps_table_data$scenario)
samps_table_data$scenario <- factor(samps_table_data$scenario, 
                                    levels(samps_table_data$scenario)[c(1,5,2,6,3,4,7)])

samps_table_data$ckmr_ssmult <- factor(samps_table_data$ckmr_ssmult, levels = c("50%","100%","150%"))

names(samps_table_data) <- c("Scenario", "Yrs", "Nsamps", "Females", "Males", "POS",
                             "POD", "HSS", "HSD", "Unrelated")
samps_table_data$Yrs <- as.character(samps_table_data$Yrs)
samps_table_data$Yrs <- ifelse(samps_table_data$Yrs == "03 Yrs", 3, 10)

samps_table_data <- samps_table_data %>% arrange(Scenario, Yrs, Nsamps)




