
## ---------------------------
##
## Script name: results_figures.R
##
## Purpose of script: generate figures for study results
##
## Author: Claudia & Eric
##
## Date Created: 2024-01-16
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

# simulated population
sim_pop <- readRDS(file.path(here(), path, "sim_pop_final.rds"))

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

layout1 <- "
ABC
ADE
"

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


# ---------------------------------------------
#
# Paper Main Text Figures
#
# ---------------------------------------------


#--------------------------------
# Figure 3: simulated reproductive value and sex-ratio-dependent skipped spawning

plot_age <- seq(1,33,0.1)
plot_age_Male <- seq(4,33,0.1)

#### fecundity

fec_base <- fec_scalar * (plot_age)^ffec_exp # females

# male fecundity at age
male_base <- ((plot_age_Male)/maxAge)^ffec_exp # base case, male fec/RV = female
male_high <- ((plot_age_Male)/maxAge)^mfec_exp_high # males steeper
male_low <- ((plot_age_Male)/maxAge)^mfec_exp_low # males flatter

male_fec_table <- data.frame(age = plot_age_Male,
                         male_base = male_base,
                         male_high = male_high,
                         male_low = male_low) %>%
  gather(key, value, -age) %>%
  separate(key, into = c("sex", "type"), sep = "_") %>%
  mutate(type = case_when(
    type == "base" ~ "mERRO = fERRO",
    type == "high" ~ "mERRO flatter",
    type == "low" ~ "mERRO steeper"
  ))

fem_fec_table <- data.frame(age = plot_age, value = fec_base/max(fec_base))  

p_fec <- ggplot(fem_fec_table) +
  geom_line(aes(x = age, y = value), linetype = "dashed", linewidth = 0.8) +
  geom_line(data = male_fec_table,
            aes(x = age, y = value, color = type)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(linetype = guide_legend(title = "mERRO OM scenario")) +
  xlab("Age") + ylab("Relative Reproductive Success")

#### sperm limitation

prop_male <- seq(0,1,0.01)
resp <- lam_fun(prop_male, 1, lfpar)

sl <- data.frame(PropMale = prop_male, Response = resp)

sl_plot <- ggplot(sl) +
  geom_line(aes(x = PropMale, y = 1-Response)) +
  theme_bw() + ylab("Probability of Skipping Spawning") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("Relative Proportion Male")

pspawn_sigmoid <- data.frame(age = ages,
                             sig = 1/(1+exp(-pspawn_steep*(ages-pspawn_infl)))
)

sig_plot <- ggplot(pspawn_sigmoid) +
  geom_line(aes(x = age, y = sig)) +
  ylab("Probability")

mRC_plot <- p_fec + sl_plot + plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A')
ggsave(here("manuscript", "figures", "Figure_3.pdf"), plot = mRC_plot, width = 7, height = 3)


#--------------------------------
# Figure 4: simulated pedigrees

sim_pop$scenario <- as.factor(sim_pop$scenario)

sim_pop2 <- sim_pop %>% 
  filter(scenario %in% c("S0","OM_S2","OM_S3", "OM_S4"),
         Year >= fyrplot) %>%
  dplyr::select(Year, scenario, Males, Females, id) %>%
  gather(data, number, Females:Males) 

sim_pop_sum <- sim_pop2 %>%
  group_by(scenario, Year, data) %>%
  summarize(median_sim = median(number),
            q25 = quantile(number, 0.25),
            q75 = quantile(number, 0.75)) %>%
  ungroup()

sim_pop_plot <- ggplot(sim_pop2) +
  geom_line(aes(x = Year, y = number, col = scenario, linetype = as.factor(id)), alpha = 0.1)  +
  geom_line(data = sim_pop_sum, 
    aes(x = Year, y = median_sim, col = scenario), linewidth = 1.5) +
  scale_linetype_manual(values = rep(1, 200)) +
  scale_color_manual(values = c("tomato3","#00BFC4", "#7CAE00", "grey30")) +
  scale_y_continuous(labels = scientific_format()) +
  facet_wrap(~data, scales = "free") +
  ylab("Abundance (Age 3+)") +
  theme(legend.position = "none")

ggsave(here("manuscript", "figures", "Figure_4.pdf"), plot = sim_pop_plot, width = 8, height = 5)


#--------------------------------
# Figure 5: ckmr pairs


ckmr_pairs_by_scenario <- ckmr_obs_samps_long %>% 
  filter(em_sc == "Base", nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, et == base_et)

# Eric explores using TeX() on changed factor levels. See: https://stackoverflow.com/questions/56518893/annotate-ggplot2-face-labels-with-latex-in-r
ckmr_pairs_by_scenario_TEX <- ckmr_pairs_by_scenario
levels(ckmr_pairs_by_scenario_TEX$pair_type) <- c(
  POD = TeX("POP$^\\neq$"),
  POS = TeX("POP$^=$"),
  HSD = TeX("HSP$^\\neq$"),
  HSS = TeX("HSP$^=$"),
  Unrelated = TeX("Unrelated")
)

main_ckmr_pairs_plot <- ggplot(ckmr_pairs_by_scenario_TEX) +
  geom_boxplot(aes(x = om_sc, y = count, fill = om_sc)) +
  facet_wrap(~pair_type, labeller = label_parsed,
             scales = "free", ncol = 3) +
  scale_fill_discrete(
    labels = om_labels
  ) +
  ylab("Number of Pairs") + xlab("") +
  theme(axis.text.x = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.70, 0.25),
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        strip.text.x = element_text(size = 16)) +
  guides(fill = guide_legend(title = "OM Scenario"))

ggsave(here("manuscript", "figures", "Figure_5.pdf"), plot = main_ckmr_pairs_plot, width = 8, height = 6)


#--------------------------------
# Figure 6: PRE by number of CKMR samples and number of sampling yrs


err_ckmr_sampling <- abund_res %>%
  filter(om_sc == "Base", em_sc == "Base") %>%
  group_by(data, et, ckmr_ssmult, nsampyrs, id) %>%
  summarize(perror = median(perror)) %>%
  ungroup()


sp1 <- err_ckmr_sampling %>% filter(et == "fixed trans & mRV") %>%
  ggplot(aes(x = ckmr_ssmult, y = perror, fill = nsampyrs)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
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
               outlier.shape = 21, outlier.stroke = 0.3) +
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
               outlier.shape = 21, outlier.stroke = 0.3) +
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
               outlier.shape = 21, outlier.stroke = 0.3) +
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
        axis.title = element_blank())

ggsave(here("manuscript", "figures", "Figure_6.pdf"), plot = ckmr_sampling_plot, width = 9, height = 9)

### stats for text:

# tmp <- err_ckmr_sampling %>%
#   filter(et != "fixed trans & mRV") %>%
#   group_by(data, et, ckmr_ssmult, nsampyrs) %>%
#     summarize(med = median(perror),
#               iqr = IQR(perror),
#               min = min(perror),
#               max = max(perror)) %>%
#     arrange(data, med) %>%
#     mutate_if(is.numeric, round, digits = 3)
# summary(tmp$med)
# summary(tmp$iqr)
#  
# err_ckmr_sampling %>%
#   filter(et != "fixed trans & mRV") %>%
#   group_by(nsampyrs, et, data, ckmr_ssmult) %>%
#   summarize(metric = IQR(perror)) %>%
#   arrange(nsampyrs, et, data, ckmr_ssmult) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   spread(nsampyrs, metric) %>%
#   mutate(diff = `10 Yrs` - `03 Yrs`) %>%
#   arrange(diff)
#  
# err_ckmr_sampling %>%
#   filter(et != "fixed trans & mRV") %>%
#   group_by(ckmr_ssmult, et) %>%
#   summarize(metric = IQR(perror)) %>%
#   arrange(ckmr_ssmult) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   spread(ckmr_ssmult, metric) %>%
#   mutate(diff50_100 = `50%` - `100%`,
#          diff100_150 = `100%` - `150%`,
#          perc50_100 = diff50_100*100/`100%`,
#          perc100_150 = diff100_150*100/`150%`) %>%
#   arrange(perc50_100)
# 
# err_ckmr_sampling %>%
#   group_by(et, ckmr_ssmult) %>%
#   summarize(metric = IQR(perror)) %>%
#   arrange(ckmr_ssmult) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   spread(ckmr_ssmult, metric) %>%
#   mutate(diff50_100 = `50%` - `100%`,
#          diff100_150 = `100%` - `150%`,
#          perc50_100 = diff50_100*100/`100%`,
#          perc100_150 = diff100_150*100/`150%`) %>%
#   arrange(perc50_100)
# 
# err_ckmr_sampling %>%
#   group_by(et, ckmr_ssmult, data) %>%
#   summarize(metric = median(perror)) %>%
#   arrange(ckmr_ssmult) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   spread(ckmr_ssmult, metric) %>%
#   mutate(diff50_100 = `50%` - `100%`,
#          diff100_150 = `100%` - `150%`,
#          perc50_100 = diff50_100*100/`100%`,
#          perc100_150 = diff100_150*100/`150%`) %>%
#   arrange(perc50_100)
# 
# err_ckmr_sampling %>%
#   group_by(et, ckmr_ssmult, data, nsampyrs) %>%
#   summarize(metric = IQR(perror)) %>%
#   ungroup() %>%
#   arrange(ckmr_ssmult) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   spread(ckmr_ssmult, metric) %>%
#   mutate(diff50_100 = `50%` - `100%`,
#          diff100_150 = `100%` - `150%`,
#          perc50_100 = diff50_100*100/`100%`,
#          perc100_150 = diff100_150*100/`150%`) %>%
#   group_by(et) %>%
#   summarize(diff50_100 = mean(diff50_100),
#             diff100_150 = mean(diff100_150),
#             perc50_100 = mean(perc50_100),
#             perc100_150 = mean(perc100_150))



#--------------------------------
# Figure 7: PRE by CKMR NLL Issues


err_emsc1_2 <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, om_sc == "Base") %>%
  group_by(em_sc, et, data, id) %>%
  summarize(perror = median(perror))

em_p1 <- 
  ggplot(err_emsc1_2 %>%
           filter(em_sc %in% "Base"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  ggtitle("Base EM") +
  ylab("Percent Relative Error") +
  coord_cartesian(ylim = c(-75,100)) +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "A")

em_p2 <- 
  ggplot(err_emsc1_2 %>%
           filter(em_sc %in% "HSPs only"), 
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
  coord_cartesian(ylim = c(-75,100)) +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "B")

em_p3 <- 
  ggplot(err_emsc1_2 %>%
           filter(em_sc %in% "est Hmt < sim Hmt"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  coord_cartesian(ylim = c(-65, 750)) +
  ggtitle(expression(bold(H[mt]^{"EM"} < H[mt]^{"OM"}))) +
  ylab("Percent Relative Error") +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "C")

em_p4 <- 
  ggplot(err_emsc1_2 %>%
           filter(em_sc %in% "sim Hmt < est Hmt"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  coord_cartesian(ylim = c(-75,100)) +
  ggtitle(expression(bold(H[mt]^{"OM"} < H[mt]^{"EM"}))) +
  ylab("Percent Relative Error") +
  guides(fill = guide_legend(title = "EM Scenario"))  +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  labs(tag = "D")

ckmr_nll_plot <- (ytitle + em_p1 + em_p2 + em_p3 + em_p4 +
                    plot_layout(
                      guides = "collect",
                      design = layout1,
                      widths = c(0.05, 0.475, 0.475)
                    ) &
                    theme(
                      legend.position = "bottom",
                      legend.justification = "center",  # <-- This centers the legend
                      axis.line.x = element_line(color = "grey20"),
                      axis.ticks.x = element_line(color = "grey20"),
                      axis.title = element_blank(),
                      axis.text.x = element_blank()
                    ) &
                    guides(fill = guide_legend(title = "EM Scenario", nrow = 2, byrow = TRUE))
                  ) +
  plot_annotation()  # Ensures the patchwork layout processes all theme elements

ggsave(here("manuscript", "figures", "Figure_7.pdf"), plot = ckmr_nll_plot, width = 9, height = 9)

# err_emsc1_2 %>% 
#   filter(em_sc %in% c("est Hmt < sim Hmt")) %>%
#   group_by(em_sc, et) %>% 
#   summarize(med = median(perror),
#             iqr = IQR(perror)) %>% 
#   arrange(et)
# 
# err_emsc1_2 %>% 
#   filter(!em_sc %in% c("HSPs only", "Base")) %>%
#   group_by(em_sc, et, data) %>% 
#   summarize(metric = median(perror)) %>% 
#   spread(em_sc, metric) %>%
#   arrange(data)
# 
# err_emsc1_2 %>% 
#   filter(!em_sc %in% c("HSPs only", "Base")) %>%
#   group_by(em_sc, et, data) %>% 
#   summarize(metric = IQR(perror)) %>% 
#   spread(em_sc, metric) %>%
#   arrange(data)


#--------------------------------
# Figure 8: Sex ratio PRE plot

sr_error <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult,
         em_sc != "HSPs only") %>%
  dplyr::select(ped_rep, nsampyrs, ckmr_ssmult,
                et, year, ckmr_seed, 
                data, sim_value, est_value, id, em_sc, om_sc) %>%
  mutate(data2 = data) %>%
  spread(data, est_value) %>%
  rename(fem_est = Females, male_est = Males) %>%
  spread(data2, sim_value) %>%
  rename(fem_sim = Females, male_sim = Males) %>%
  group_by(ped_rep, nsampyrs, ckmr_ssmult,
           et, year, ckmr_seed,
           id, em_sc, om_sc) %>%
  summarize(fem_est = sum(fem_est, na.rm = TRUE),
            male_est = sum(male_est, na.rm = TRUE),
            fem_sim = sum(fem_sim, na.rm = TRUE),
            male_sim = sum(male_sim, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(sim_sr = male_sim/(male_sim+fem_sim)*100,
         est_sr = male_est/(male_est+fem_est)*100,
         prerror = get_perror(sim = sim_sr, est = est_sr))

sr_error_sum <- sr_error %>%
  group_by(om_sc, em_sc, et, id) %>%
  summarize(prerror = median(prerror))

# sr_error_sum <- sr_error %>%
#   droplevels() %>%
#   group_by(ckmr_ssmult,scenario,nsampyrs,
#            year
#   ) %>%
#   mutate(N = n()) %>%
#   ungroup() %>%
#   group_by(ckmr_ssmult,scenario,nsampyrs) %>%
#   mutate(propN = N/max(N)) %>%
#   ungroup() %>%
#   filter(propN > prop2plot)



# sr_error_sum %>%
#   filter(em_sc == "est Hmt < sim Hmt") %>%
#   group_by(et) %>%
#   summarize(med = median(prerror),
#             iqr = IQR(prerror),
#             min = min(prerror),
#             max = max(prerror)) %>%
#   mutate_if(is.numeric, round, digits = 3) %>%
#   arrange(et)

sryll <- -95
sryul <- 150

layout3 <- "
ABCD
AEFG
"


sr1 <- ggplot(sr_error_sum %>% filter(om_sc == "Base", em_sc == "Base"), 
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(sryll,sryul)) +
  ggtitle("Base") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "A")

sr2 <- ggplot(sr_error_sum %>% filter(om_sc == "skipped spawning"), 
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(sryll,sryul)) +
  ggtitle("Skipped Spawning") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        #axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "B")

sr3 <- ggplot(sr_error_sum %>% filter(em_sc == "est Hmt < sim Hmt"), 
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(sryll,sryul)) +
  ggtitle(expression(bold(H[mt]^{"EM"} < H[mt]^{"OM"}))) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "C")

sr4 <- ggplot(sr_error_sum %>% filter(em_sc == "sim Hmt < est Hmt"), 
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(-5, 310)) +
  ggtitle(expression(bold(H[mt]^{"OM"} < H[mt]^{"EM"}))) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        #axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "D")

sr5 <- ggplot(sr_error_sum %>% filter(om_sc == "mRV flatter"), 
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(sryll,sryul)) +
  ggtitle(expression(bold(RV[M]^{"OM"} ~ "flatter than" ~ RV[F]^{"OM"}))) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "E")

sr6 <- ggplot(sr_error_sum %>% filter(om_sc == "mRV steeper"), # 533 omitted by coord_cart 200
              aes(x = et, y = prerror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  geom_boxplot(aes(fill = et), outlier.size = 0.6, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3,
               linewidth = 0.2, fatten = 5) +
  coord_cartesian(ylim = c(sryll,sryul)) +
  ggtitle(expression(bold(RV[M]^{"OM"} ~ "steeper than" ~ RV[F]^{"OM"}))) +
  scale_fill_discrete(
    labels = em_labels
  ) +
  guides(fill = guide_legend(title = "Parameter Estimation")) +
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        #axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(tag = "F")

sr_plot <- ytitle + sr1 + sr2 + sr3 + sr4 + sr5 + sr6 +
  plot_layout(
    axis = "collect",
    guides = "collect",
    design = layout3,
    widths = c(0.05, 0.315, 0.315, 0.32)
  ) &
  theme(
    axis.line.x = element_line(color = "grey20"),
    axis.ticks.x = element_line(color = "grey20"),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  ) &
  guides(fill = guide_legend(title = "EM Scenario", nrow = 2, byrow = TRUE))

sr_plot_final <- sr_plot +
  plot_annotation() &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )


ggsave(here("manuscript", "figures", "Figure_8.pdf"), plot = sr_plot_final, width = 8, height = 8)


#--------------------------------
# Figure 9: Estimated function shapes


est_tf <- sim_pars %>% 
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, em_sc != "HSPs only",
         par_name == "TFsd") %>% 
  dplyr::select(id, om_sc, et, value, sim_value, em_sc) %>%
  rename(est_par = value, sim_par = sim_value) %>%
  mutate(sc3 = case_when(
    em_sc == "est Hmt < sim Hmt" ~ "eHmt < sHmt",
    em_sc == "sim Hmt < est Hmt" ~ "sHmt < eHmt",
    om_sc == "Base" ~ "Base",
    TRUE ~ om_sc
  ))
est_tf$id2 <- 1:nrow(est_tf)
est_tf$sc3 <- factor(est_tf$sc3, 
                     levels = c("Base",
                                "skipped spawning",
                                "eHmt < sHmt",
                                "sHmt < eHmt",
                                "mRV flatter",
                                "mRV steeper"))


est_tf <- est_tf %>%
  filter(et %in% c("est trans & fixed mRV", "est trans & mRV")) %>%
  mutate(et = ifelse(et == "est trans & fixed mRV", "mRV fixed", "mRV estimated")) %>%
  rowwise() %>%
  mutate(est_tfun = list(pnorm(q = 1:33, mean = t1, sd = est_par))) %>%
  unnest_wider(est_tfun, names_sep = "_") %>%
  pivot_longer(cols = starts_with("est_tfun_"), names_to = "Age", values_to = "est_value")
est_tf$Age <- as.numeric(gsub("est_tfun_", "", est_tf$Age))
est_tf$sim_value <- rep(pnorm(q=1:33, mean=t1, sd=t2), nrow(est_tf)/33)

# Plotting
st_plot <- ggplot(est_tf) +
  geom_line(aes(x = Age, y = est_value, group = id2, color = sc3), alpha = 0.1) +
  geom_line(data = est_tf %>%
              distinct(et, sc3, Age, sim_value),
            aes(x = Age, y = sim_value), linetype = "dashed", linewidth = 0.8) +
  geom_line(data = est_tf %>%
              group_by(et, sc3, Age) %>% summarize(median_val = median(est_value)),
            aes(x = Age, y = median_val), linewidth = 0.8) +
  facet_grid(sc3~ et) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "none",
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ylab("Probability") + ggtitle("Estimated Sex Transition Fct")


est_mf <- sim_pars %>% 
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, ckmr_ssmult == base_ckmr_ssmult,
         em_sc != "HSPs only",
         par_name == "mRV_exp") %>% 
  dplyr::select(id, om_sc, et, value, sim_value, em_sc) %>%
  rename(est_par = value, sim_par = sim_value) %>%
  mutate(sc3 = case_when(
    em_sc == "est Hmt < sim Hmt" ~ "eHmt < sHmt",
    em_sc == "sim Hmt < est Hmt" ~ "sHmt < eHmt",
    om_sc == "Base" ~ "Base",
    TRUE ~ om_sc
  )) 
est_mf$id2 <- 1:nrow(est_mf)
est_mf$sc3 <- factor(est_mf$sc3, 
                     levels = c("Base",
                                "skipped spawning",
                                "eHmt < sHmt",
                                "sHmt < eHmt",
                                "mRV flatter",
                                "mRV steeper"))

est_mf1 <- est_mf %>%
  filter(et %in% c("fixed trans & est mRV", "est trans & mRV")) %>%
  mutate(et = ifelse(et == "fixed trans & est mRV", "transition fixed", "transition estimated")) %>%
  rowwise() %>%
  mutate(est_mfun = list(((1:33)/33)^est_par)) %>%
  unnest_wider(est_mfun, names_sep = "_") %>%
  pivot_longer(cols = starts_with("est_mfun_"), names_to = "Age", values_to = "est_value")
est_mf1$Age <- as.numeric(gsub("est_mfun_", "", est_mf1$Age))

est_mf2 <- est_mf %>%
  filter(et %in% c("fixed trans & est mRV", "est trans & mRV")) %>%
  mutate(et = ifelse(et == "fixed trans & est mRV", "transition fixed", "transition estimated")) %>%
  rowwise() %>%
  mutate(est_mfun = list(((1:33)/33)^sim_par)) %>%
  unnest_wider(est_mfun, names_sep = "_") %>%
  pivot_longer(cols = starts_with("est_mfun_"), names_to = "Age", values_to = "sim_value")
est_mf2$Age <- as.numeric(gsub("est_mfun_", "", est_mf2$Age))

est_mf <- est_mf1 %>%
  left_join(est_mf2)


me_plot <- ggplot(est_mf) +
  geom_line(aes(x = Age, y = est_value, group = id2, color = sc3), alpha = 0.1) +
  geom_line(data = est_mf %>%
              distinct(et, sc3, Age, sim_value),
            aes(x = Age, y = sim_value), linetype = "dashed", linewidth = 0.8) +
  geom_line(data = est_mf %>%
              group_by(Age, sc3, et) %>% summarize(median_val = median(est_value)),
            aes(x = Age, y = median_val), linewidth = 0.8) +
  facet_grid(sc3 ~ et) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ylab("Reproductive Value") + ggtitle("Estimated Male RV Fct")


est_fun_plot <- st_plot + me_plot + plot_layout(axis_titles = "collect_x", axes = "collect") &
  theme(axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"))

ggsave(here("manuscript", "figures", "Figure_9.pdf"), plot = est_fun_plot, width = 11, height = 12)


#--------------------------------
# Figure 10: PRE by male reproductive value scenarios

mr_res <- abund_res %>%
  filter(nsampyrs == base_ckmr_nsampyrs, sample_by_sex == samp_sex_base, 
         ckmr_ssmult == base_ckmr_ssmult, em_sc == "Base") %>%
  group_by(om_sc, et, data, id) %>%
  summarize(perror = median(perror)) %>%
  ungroup()

mrole_y_upper <- 100

om_p1 <-
  ggplot(mr_res %>% filter(om_sc == "Base"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_x_discrete(labels = function(x) str_replace_all(x, "&", "&\n")) +
  ylab("Percent Relative Error") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-75,mrole_y_upper)) +
  guides(fill = guide_legend(title = "EM Scenario")) +
  ggtitle("Base OM") +
  labs(tag = "A")

om_p2 <-
  ggplot(mr_res %>% filter(om_sc == "skipped spawning"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_x_discrete(labels = function(x) str_replace_all(x, "&", "&\n")) +
  ylab("Percent Relative Error") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-75,mrole_y_upper)) + # omits 26
  guides(fill = guide_legend(title = "EM Scenario")) +
  ggtitle("Skipped spawning") +
  labs(tag = "B")

om_p3 <-
  ggplot(mr_res %>% filter(om_sc == "mRV flatter"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) +
  scale_x_discrete(labels = function(x) str_replace_all(x, "&", "&\n")) +
  ylab("Percent Relative Error") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-75,mrole_y_upper)) +
  guides(fill = guide_legend(title = "EM Scenario")) +
  ggtitle(expression(bold(RV[M]^{"OM"} ~ "flatter than" ~ RV[F]^{"OM"}))) +
  labs(tag = "C")

om_p4 <-
  ggplot(mr_res %>% filter(om_sc == "mRV steeper"), 
         aes(x = et, y = perror)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_boxplot(aes(fill = et),
               outlier.size = 0.8, outlier.fill = NULL, 
               outlier.shape = 21, outlier.stroke = 0.3) +
  facet_wrap(~data) + 
  scale_x_discrete(labels = function(x) str_replace_all(x, "&", "&\n")) +
  ylab("Percent Relative Error") +
  scale_fill_discrete(
    labels = em_labels
  ) +
  theme(axis.text = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(-75,mrole_y_upper)) +
  guides(fill = guide_legend(title = "EM Scenario")) +
  ggtitle(expression(bold(RV[M]^{"OM"} ~ "steeper than" ~ RV[F]^{"OM"}))) +
  labs(tag = "D")

male_role_plot <- ytitle + om_p1 + om_p2 + om_p3 + om_p4 +
  #plot_annotation(tag_levels = 'A') +
  plot_layout(#axis_titles = "collect",
    guides = "collect",
    design = layout1,
    widths = c(0.05,0.475,0.475)) &
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "grey20"),
        axis.ticks.x = element_line(color = "grey20"),
        axis.title = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        #strip.text.y = element_blank(),
        axis.text.x = element_blank()) &
  guides(fill = guide_legend(title = "EM Scenario", nrow = 2,
                             byrow = T))

male_role_plot_final <- male_role_plot +
  plot_annotation() &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

ggsave(here("manuscript", "figures", "Figure_10.pdf"), plot = male_role_plot_final, width = 9, height = 9)


# mr_res %>%
#   filter(om_sc == "Base", data == "Females") %>%
#   summarize(med = median(perror))

