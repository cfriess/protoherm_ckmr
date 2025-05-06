localrules: estimate


# Just a note on running this.  For making the pedigrees, spread it across the
# himem and medmem machines:
#snakemake -np --profile hpcc-profiles/slurm/jdb/sedna --use-envmodules  --rerun-triggers mtime --until make_pedigree --keep-going

# to do the estimate jobs, get a medmem machine checked out and then do:
# snakemake -np --cores 20 --use-envmodules  --rerun-triggers mtime --keep-going



# Here is a function with reasonable defaults for setting up different
# simulation targets
def sim_targets(
  pop_mult=256,
  mating="random",    #["random","cmd"],
  ssf=0.01,           #[0.01, 0.8],                                 
  eofsr=0,            #[0, 1, 2],                           
  nsampyrs=[3,10],
  ped_rep=[f"{i:03}" for i in range(1, 40+1)],
  #ped_rep=[f"{i:03}" for i in range(1, 4+1)],
  ckmr_ssmult=[0.5,1,1.5],
  rind_sd=0.01, 
  prop_sampled=1, 
  fec_known=1,   #[1,2,3]
	beta=1, 
  comp_bias="FALSE",
  hzt=1,              #[0.9,1],
  hze=1,              #[0.9,1],
  et=[0,1],
  sample_by_sex=1,
  ckmr_seed=[f"{i:03}" for i in range(1, 10+1)],
  ri_seed=[f"{i:03}" for i in range(1, 2+1)],
  nll_type=1,         #[0,1],
  cv_seed=[f"{i:03}" for i in range(1, 2+1)]):
    return expand("results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/analysis_objs.rds",
			pop_mult=pop_mult,
			mating=mating,
			ssf=ssf,                                 # could be good to play with more
			eofsr=eofsr,                           # TRUE is closer to the estimation model 
			nsampyrs=nsampyrs,
			ped_rep=ped_rep,
			ckmr_ssmult=ckmr_ssmult,
			rind_sd=rind_sd,
			prop_sampled=prop_sampled,
			fec_known=fec_known,
			beta=beta,
			comp_bias=comp_bias,
			hzt=hzt,
			hze=hze,
			et=et,
			sample_by_sex=sample_by_sex,
			ckmr_seed=ckmr_seed,
			ri_seed=ri_seed,
			cv_seed=cv_seed,
			nll_type=nll_type
			)


### ACTUAL TARGETS SIMULATED:

### For Exploratory runs (life history factors):
# 40 peds, 10 samples from each
# we don't need ri_seed and cv_seeds with the default settings 
# all other values = defaults
### Exploration runs for fishery data uncertainty and bias (no need to remake pedigrees)
### See MS at 544b6340d55 for terminology
s0=sim_targets(ri_seed=1, cv_seed=1)  # base - for 3 levels of nsampyrs, 2 levels of ckmr_ssmult and 2 levels of et
s1=sim_targets(mating="cmd", ssf=0.8, ri_seed=1, cv_seed=1) #  harem spawning and high fem spawning group fidelity
s2=sim_targets(eofsr=3, ri_seed=1, cv_seed=1)  # skipped spawning, younger females have higher prob of skipping
s3=sim_targets(fec_known=2,ri_seed=1, cv_seed=1)  # male rel repro output higher than EM
s4=sim_targets(fec_known=3, ri_seed=1, cv_seed=1)  # male rel repro output lower than EM

e1=sim_targets(ri_seed=1, cv_seed=1, hze=0.9)  # hze < hzt
e2=sim_targets(ri_seed=1, cv_seed=1, hzt=0.9)  # hzt < hze
e3=sim_targets(ri_seed=1, cv_seed=1, nll_type=0)  # HSPs only in likelihood
e4=sim_targets(ri_seed=1, cv_seed=1, et=2) # fixed transition, estimate maleFec
e5=sim_targets(ri_seed=1, cv_seed=1, et=3) # estimate both transition and maleFec
e6=sim_targets(ri_seed=1, cv_seed=1, fec_known=2, et=2) # fixed transition, estimate maleFec
e7=sim_targets(ri_seed=1, cv_seed=1, fec_known=2, et=3) # estimate both transition and maleFec
e8=sim_targets(ri_seed=1, cv_seed=1, fec_known=3, et=2) # fixed transition, estimate maleFec
e9=sim_targets(ri_seed=1, cv_seed=1, fec_known=3, et=3) # estimate both transition and maleFec

e10=sim_targets(ri_seed=1, cv_seed=1, hze=0.9, et=2) # hze < hzt, estimate maleFec
e11=sim_targets(ri_seed=1, cv_seed=1, hzt=0.9, et=2) # hzt < hze, estimate maleFec
e12=sim_targets(ri_seed=1, cv_seed=1, hze=0.9, et=3) # hze < hzt, estimate transition and maleFec
e13=sim_targets(ri_seed=1, cv_seed=1, hzt=0.9, et=3) # hzt < hze, estimate transition and maleFec
e14=sim_targets(ri_seed=1, cv_seed=1, nll_type=0, et=2)  # HSPs only, estimate maleFec
e15=sim_targets(ri_seed=1, cv_seed=1, nll_type=0, et=3)  # HSPs only, estimate transition and maleFec
e16=sim_targets(ri_seed=1, cv_seed=1, et=2, eofsr=3) # fixed trans, est maleFec, skipped spawning
e17=sim_targets(ri_seed=1, cv_seed=1, et=3, eofsr=3) # est transition and maleFec, skipped spawning

e18=sim_targets(ri_seed=1, cv_seed=1, eofsr=3, nll_type=[0,1], et = [0,1,2,3])  # HSPs only in likelihood for skipped spawning

e19=sim_targets(ri_seed=1, cv_seed=1, et = [0,1,2,3], sample_by_sex = [2,3], eofsr=[0,3], nsampyrs = 3, ckmr_ssmult = 1) # different options for # of males sampled

ALL_RUNS = s0 + s1 + s2 + s3 + s4 + e1 + e2 + e3 + e4 + e5 + e6 + e7 + e8 + e9 + e10 + e11 + e12 + e13 + e14 + e15 + e16 + e17 + e18 + e19


# Here is all the stuff we want to make then
rule all:
	input:
		ALL_RUNS

rule make_pedigree:
	params:
		pop_mult="{pop_mult}",
		mating="{mating}",
		ssf="{ssf}",
		eofsr="{eofsr}",
		fec_known="{fec_known}"
	envmodules:
		"R/4.0.3"
	log: 
		log="results/logs/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/ped_rep-{ped_rep}/ped_rlog.txt",
		snake_obj="results/logs/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/ped_rep-{ped_rep}/make_pedigree_snake_obj.txt"
	benchmark: 
		"results/benchmarks/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/ped_rep-{ped_rep}/make_pedigree.bmk"
	output:
		rdata="results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/ped_rep-{ped_rep}/ped.Rdata"
	script:
		"code/main_sim/ckmr_sim.R"


rule estimate:
	input:
		in_list="results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/ped_rep-{ped_rep}/ped.Rdata"
	params:
	  rind_sd="{rind_sd}",
	  prop_sampled="{prop_sampled}",
	  nsampyrs="{nsampyrs}",
	  ckmr_ssmult="{ckmr_ssmult}",
	  fec_known="{fec_known}",
	  beta="{beta}",
	  comp_bias="{comp_bias}",
	  hzt="{hzt}",
	  hze="{hze}",
	  et="{et}",
	  sample_by_sex="{sample_by_sex}",
	  ckmr_seed="{ckmr_seed}",
	  ri_seed="{ri_seed}",
	  cv_seed="{cv_seed}",
	  nll_type="{nll_type}",
	output: 
		rds="results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/analysis_objs.rds",
		key_obs="results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/key_output_objs.rds",
		admb_arena=directory("results/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/admb_arena")
	envmodules:
		"R/4.0.3"
	log: 
		log="results/logs/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/estimate_rlog.txt",
		snake_obj="results/logs/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/snake_obj.rds"
	benchmark: 
		"results/benchmarks/POP-pop_mult-{pop_mult}/MATING-mating-{mating}-ssf-{ssf}-eofsr-{eofsr}-fec_known-{fec_known}/NSAMPYRS-{nsampyrs}/ped_rep-{ped_rep}/ckmr_ssmult-{ckmr_ssmult}/EST-rind_sd-{rind_sd}-prop_sampled-{prop_sampled}-beta-{beta}-comp_bias-{comp_bias}/EST2-mtHZ-{hzt}-mtHZ_est-{hze}-et-{et}/sample_by_sex-{sample_by_sex}---ckmr_seed-{ckmr_seed}---ri_seed-{ri_seed}---cv_seed-{cv_seed}---nll_type-{nll_type}/estimate.txt",
	script:
		"code/main_sim/estimate-SMK.R"
