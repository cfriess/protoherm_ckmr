
#### table 1 in paper

library(kableExtra)

col1 <- c(
  '$t$',
  '$a$',
  '$A$',
  '$i$,$j$',
  '$c$',
  '$d$',
  'HSP',
  'POP',
  '$a^*$',
  '$\\mathrm{ERRO}_{a,t,g}$',
  '$g_i$ or $g_{i,s_i}$',
  '$N_a$',
  '$N_a^*$',
  '$\\mathrm{Mat}_a$',
  '$\\mathrm{RV}_{a,F}^\\mathrm{OM}$, $\\mathrm{RV}_{a,F}^\\mathrm{EM}$',
  '$\\mathrm{RV}_{a,M}^\\mathrm{OM}$, $\\mathrm{RV}_{a,M}^\\mathrm{EM}$',
  '$S_a$',
  '$\\mathrm{HSP}_{F\\rightarrow F}$, $\\mathrm{HSP}_{F\\rightarrow M}$, $\\mathrm{HSP}_{M\\rightarrow M}$', 
  '$P_a^{\\mathrm{F}\\rightarrow \\mathrm{M}}$,',
  'mtDNA',
  '$\\mathrm{mt}_{a,b}$',
  '$H_\\mathrm{mt}$',
  '$H_\\mathrm{mt}^\\mathrm{OM}$, $H_\\mathrm{mt}^\\mathrm{EM}$',
  '$s$',
  '$\\mathrm{POP}_{F\\rightarrow F}$, $\\mathrm{POP}_{F\\rightarrow M}$, $\\mathrm{POP}_{M\\rightarrow M}$', 
  '$\\mathrm{ETRO}_{t,g}$',
  '$\\K{HSP}{c_i,c_j}{g_i,s_i}{a,b}$',
  '$\\K{POP}{c_i,c_j}{g_i,s_i}{a,b}$',
  '$\\K{Not}{c_i,c_j}{g_i,s_i}{a,b}$',
  '$M_a$',
  '$v_{a,f}$', 
  '$F_{t,f}$',
  '$n_\\mathrm{base}$',
  '$\\rbar$',
  '$\\mathrm{STF}_{SD}$',
  '$\\mathrm{PRE}$',
  '$\\mathrm{IQR}$'
)

col2 <- c(
  "Time in years",
  "Age in years",
  "Maximum age of a fish in years",
  "Subscript denoting the first- and second-born, respectively, in a sampled pair of fish",
  "Birth year of a sampled fish",
  "The age difference between $i$ and $j$, $c_j-c_i$",
  "Half-sibling pair", 
  "Parent-offspring pair",
  "For HSPs, the shared parent's (unknown) age at the time $j$ was born",
  "Expected relative reproductive output of sex-$g$ fish in age class $a$ at time $t$, defined as the total reproductive output of that age class divided by the total reproductive output of all ages",
  "Sex of sampled fish $i$ at the time of its sampling",
  'Numbers at Age',
  'Relative Numbers at Age, meaning $N_a/\\sum{N_a}$',
  'Maturity at age',
  'Female reproductive value (at age) used in the OM and EM, respectively (for the females, $\\mathrm{RV}_F^\\mathrm{OM}$ was always equal to $\\mathrm{RV}_F^\\mathrm{EM}$)',
  'Male reproductive value (at age) simulated in the OM and used in the EM, respectively',
  "Survival at age",
  "An HSP where the shared parent was female at the time of $i$ and $j$'s births, female at $i$'s birth and male at $j$'s, or male at both $i$ and $j$'s birth, respectively.",
  "Probability that a fish that starts age $a$ as a female will switch sex to become male before reproducing",
  "Mitochondrial DNA",
  "The observation that the first-born carries a mtDNA haplotype denoted $a$ and the second-born carries a haplotype denoted $b$",
  'mitochondrial haplotype diversity---the probability that two mtDNA sequences randomly sampled from the population are different. 
  An $H_\\mathrm{mt}$ of 1 would indicate a 100\\% certainty that two HSPs with the same mtDNA shared a mother, while any value <1 means there is
a chance (which equals 1 - $H_\\mathrm{mt}$) that their mtDNA sequence is the same even though they did not
share a mother. For POPs, $H_\\mathrm{mt}$ = 1 means certainty that a parent which was male at the time of sampling was female when it produced
the offspring, if mtDNA matches',
  'Simulated $H_\\mathrm{mt}$ and $H_\\mathrm{mt}$ assumed in the estimation model, respectively',
  "Sampling year",
  "A parent offspring pair in which the parent was female when giving birth and when sampled, female when giving birth and male when sampled, or male when fertilizing the offspring and when sampled, respectively",
  "Expected total reproductive output of sex-$g$ fish at time $t$",
  "Number of half-sibling pairs identified with the first-born from time $c_i$ sampled as sex $g_i$ at time $s_i$ carrying mtDNA sequence $a$ and with the second-born from time $c_j$ carrying mtDNA sequence $b$",
  "Number of parent-offspring pairs identified with the parent born at time $c_i$ sampled as sex $g_i$ at time $s_i$ carrying mtDNA sequence $a$ and with the offspring born at time $c_j$ carrying mtDNA sequence $b$",
  "Number of unrelated sample pairs with the first-born from time $c_i$ sampled as sex $g_i$ at time $s_i$ carrying mtDNA sequence $a$ and with the second-born from time $c_j$ carrying mtDNA sequence $b$",
  'Natural mortality at age',
  'Vulnerability at age to fleet $f$',
  'Fishing mortality in year $t$ due to fleet $f$',
  'Base number of individuals sampled for CKMR, defined as 10$\\sqrt{N}$ where $N$ is mean (age 2+) abundance over the CKMR sampling period',
  'Mean Recruitment over the assessment time period, an estimated parameter',
  'Standard deviation of the $P_a^{\\mathrm{F}\\rightarrow \\mathrm{M}}$ function, an estimated parameter in some runs',
  'Percent relative error, calculated as $100(\\mathrm{EV}-\\mathrm{SV})/\\mathrm{SV}$, where $\\mathrm{EV} = \\mathrm{simulated value}$ and $\\mathrm{SV} = \\mathrm{estimated value}$',
  'Interquartile range'
)



table1 <- data.frame(Notation = col1, 
                     Definition = col2 
) 

