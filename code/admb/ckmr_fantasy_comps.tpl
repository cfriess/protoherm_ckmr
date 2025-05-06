
DATA_SECTION
	//!!CLASS ofstream MCMCreport("MCMCreport.csv",ios::trunc);
	init_adstring datafile;
	init_adstring ctlfile;
	!!ad_comm::change_datafile_name(datafile);
	init_int fyr; //first year
	init_int lyr; //last year
	init_int fage_adult; //first adult yr for females
	init_int fage; //first age group
	init_int lage; //first age group
	vector ages(fage,lage); //age vector
	!!ages.fill_seqadd(fage,1); //create a sequence, fill vector starting with fage and incrementing by 1
	init_int fckmr; //start index for ckmr samples, 1
	init_int nckmr; //number of c1 & c2 ckmr pairs
	init_matrix Strue(fyr-lage+1,lyr,fage,lage); //com proportions at age
	init_vector Af(fage,lage); //maturity at age
	init_vector Fec(fage,lage); //female fecundity at age
	//init_vector FecM(fage,lage); //male fecundity at age
	//init_vector Pfm(fage,lage); //prob of transitioning
	//init_vector Imp(fage,lage); //initial male props
	init_vector Rt(fyr-lage+1,lyr); //recruitment index
	init_number rsd; //recruitment index sd
	init_matrix comps(fyr,lyr,fage,lage); //com proportions at age
	//!!cout<<comps<<endl; 
	//!!ad_exit(1);
	init_vector vulr(fage,lage); //recreational fleet vulnerability at age
  init_vector vulc(fage,lage); //commercial fleet vulnerability at age
  init_vector Mage(fage,lage); //natural mortality at age;
	init_vector c1(fckmr,nckmr); //first-born birth year
	init_vector c2(fckmr,nckmr); //second born birth year
	init_vector s1(fckmr,nckmr); //first born sampling year
	init_vector g1(fckmr,nckmr); //first born sex
	init_vector PO(fckmr,nckmr); //number of PO pairs
	init_vector MO(fckmr,nckmr); //number of MO pairs
	init_vector HSD(fckmr,nckmr); //number of HSD pairs
	init_vector HSS(fckmr,nckmr); //number of HSS pairs
	init_vector Not(fckmr,nckmr); //number of unrelated pairs
	init_number mtHZ_est; // estimated heterozygosity
	init_int et; // binary variable controlling whether or not transition parameters are estimated
	init_int ef; // binary variable controlling whether or not male contribution is estimated
	init_number comps_weight;
	init_number index_weight;
	init_number ckmr_weight;
	init_int nll_type;
	init_int eof;
	int iter; // for mcmc output file
	!!iter=0;

	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	//load in inital parmaters NOTE that you could use a pin file see ADMB manual.
	// the pin file is useful for starting overdispersed MCMC runs 
	!!ad_comm::change_datafile_name(ctlfile);

	init_number irbar; // inital guess for average recuitment over the timeseries
	init_number ifbarr; // inital guess for average rec fishing mortality rate
	init_number ifbarc; // inital guess for average com fishing mortality rate
	init_number itinfl; // transition prob at age inflection point
	init_number itsd; // transition prob at age sd
	init_number imfec_exp; // male relative repro output exponent
	init_int eofc;

	LOCAL_CALCS
		if(eofc!=999)
		{
			cout<<"Error reading control.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS


PARAMETER_SECTION
	init_bounded_number log_rbar(6,17,1); //average recruitment over the time series
	init_bounded_number log_fbarc(-4.6,0.5,2); //average fishing mortality rate over the time series
	init_bounded_number log_fbarr(-4.6,0.5,2); //average fishing mortality rate over the time series
	init_bounded_number log_tinfl(1.5,3,-1); //transition prob at age inflection point
	init_bounded_number log_tsd(0.9,2.4,et); //transition prob at age sd
	init_bounded_number log_mfec_exp(-4,2,ef) //rel male repro output power funct exponent
	init_bounded_vector wt(fyr-lage,lyr,-10.,10.,3); //Recruitment deviation from r_bar for all inital ages and all recruitments
	init_bounded_vector log_ftc_dev(fyr-lage+1,lyr,-10.,10.,4); // Com Fishing mortality deviations from F_bar
	init_bounded_vector log_ftr_dev(fyr-lage+1,lyr,-10.,10.,4); // RecFishing mortality deviations from F_bar
	objective_function_value nll; //Total negative log likelihood
	LOCAL_CALCS
		log_rbar=log(irbar);
		log_fbarc=log(ifbarc);
		log_fbarr=log(ifbarr);
	  log_tinfl=log(itinfl);
		log_tsd=log(itsd);
		log_mfec_exp=log(imfec_exp);
	END_CALCS
	// the following are for backtransforming parameters
	number fbarc;
	number fbarr;
	number tinfl;
	number tsd;
	number mfec_exp;
	number q; //y=qVB from exp(Z_bar)

	vector lxf(fage,lage); // fished survivorship at age
	vector ftc(fyr-lage+1,lyr); // commercial fishing mortality rates
	vector ftr(fyr-lage+1,lyr); // recreational fishing mortality rates
	vector ft(fyr-lage+1,lyr); //fishing mortality rates
	matrix Z(fyr-lage+1,lyr,fage,lage); //Total mortality at age M+F*vul
	matrix S(fyr-lage+1,lyr,fage,lage); //Total survival
	matrix Frec(fyr-lage+1,lyr,fage,lage); //rec F at age and t
	matrix Fcom(fyr-lage+1,lyr,fage,lage); //com F at age and t
	matrix Nfta(fyr-lage+1,lyr+1,fage,lage);// FemaleNumbers at age
	matrix Nmta(fyr-lage+1,lyr+1,fage,lage);// Male Numbers at age
	matrix Nta(fyr-lage+1,lyr+1,fage,lage);// Combined Numbers at age
	matrix Nt(fyr,lyr+1,fage,lage);//Numbers at age
	sdreport_vector Nft(fyr,lyr+1); // summed female abundance
	sdreport_vector Nmt(fyr,lyr+1); // summed male abundance
	sdreport_vector ANft(fyr,lyr+1); // summed adult female
	sdreport_vector ANt(fyr,lyr+1); // summed total adult adults
	matrix Rf(fyr,lyr,fage,lage);// predicted female proportions at age
	matrix Rm(fyr,lyr,fage,lage);// predicted male proportions at age
	vector Pfm(fage,lage); // probability of transitioning to male
	vector Imp(fage,lage); // initial proportion male at age
	vector FecM(fage,lage); // relative male repro output
	vector NF(fage,lage); // for calculating initial proportion male, number of females
	vector NM(fage,lage); // for calculating initial proportion male, number of males
	matrix pahat(fyr,lyr,fage,lage);//predicted comps
	number Pffcc;
	number Pfmcc;
	number Pmmcc;
	number Ppopffcc;
	number Ppopfmcc;
	number Ppopmmcc
	number PHSS;
	number PHSD;
	number PPHSS;
	number PPHSD;
	vector Ef_denom(fyr,lyr);
  vector Em_denom(fyr,lyr);
  vector Pffcc_vec(fckmr,nckmr);
  vector Pfmcc_vec(fckmr,nckmr);
  vector Pmmcc_vec(fckmr,nckmr);
  vector PHSS_vec(fckmr,nckmr);
  vector PHSD_vec(fckmr,nckmr);
  vector Ppopffcc_vec(fckmr,nckmr);
  vector Ppopfmcc_vec(fckmr,nckmr);
  vector Ppopmmcc_vec(fckmr,nckmr);
  vector PPHSD_vec(fckmr,nckmr);
  vector PPHSS_vec(fckmr,nckmr);
  sdreport_vector rt_hat(fyr-lage+1,lyr+1); // predicted recruitment
  vector rt_resid(fyr-lage+1,lyr);// recruitment index residuals for likelihood

PRELIMINARY_CALCS_SECTION

PROCEDURE_SECTION
	initialization();
	statedynamics();
	fisheries_obs_model();
	ckmr_obs_model();
	objective_function();

	if(mceval_phase())
	{ 
		mcmc_output();
	}

FUNCTION initialization
  // transform parameters
  fbarc=mfexp(log_fbarc);
  fbarr=mfexp(log_fbarr);
  tinfl=mfexp(log_tinfl);
	tsd=mfexp(log_tsd);
  mfec_exp=mfexp(log_mfec_exp);
  
  //fishing mortalites on fully vunerable
	ftc=mfexp(log_fbarc+log_ftc_dev);	// adding lognormal errors
	ftr=mfexp(log_fbarr+log_ftr_dev);	
	Fcom=outer_prod(ftc,vulc); // commercial fishing mortality rates by age
	Frec=outer_prod(ftr,vulr); // recreational fishing mortality rates by age
	for(int i=fyr-lage+1;i<=lyr;i++)
	{
		Z(i)=Mage;
	}
	for(int i=fyr-lage+1;i<=lyr;i++)
	{
		Z(i)+=Fcom(i)+Frec(i);
		S(i)=mfexp(-Z(i));
	}
	
	// fished incidence function
	lxf(1)=1.;
	for(int i=(fage+1);i<=lage;i++)
	{
		lxf(i)=lxf(i-1)*S(fyr-lage+1,i-1);
	}
	lxf(lage)/=(1.-1.*S(fyr-lage+1,lage));//add the plus group
	
	// get initial proportion male
	// and male repro output
	for(int a=fage;a<=lage;a++)
	{
		NF(a)=1.0;
		NM(a)=0.0;
		Pfm(a)=cumd_norm((ages(a)-tinfl)/tsd);
		FecM(a) = pow(static_cast<double>(a) / lage, mfec_exp);
	}
	//cout<<FecM<<endl;
	//ad_exit(1);
  for(int i=2;i<=lage;i++){ // 
  	NM(i)=NM(i-1)+NF(i-1)*Pfm(i-1);
  	NF(i)=NF(i-1)*(1-Pfm(i-1));
  	Imp(i)=NM(i)/(NM(i)+NF(i));
  }
  for(int a=1;a<=3; a++) {
    Pfm(a)=0.0; // no 3-year old males
    Imp(a)=0.0; // no 3-year old males
  }
	
	// initialize numbers at age in yr fyr-lage+1
	Nfta(fyr-lage+1,fage)=mfexp(log_rbar+wt(fyr-lage));//initalize numbers at age 1 in the first year as rbar and deviation
	Nmta(fyr-lage+1,fage)=0.0; // no 1 yr old males
	rt_hat(fyr-lage+1)=Nfta(fyr-lage+1,fage); // first yr recruitment
	for(int j=(fage+1);j<=lage;j++)
	{
	  Nfta(fyr-lage+1,j)=mfexp(log_rbar+wt(fyr-lage))*lxf(j);//first yr numbers at age as rbar and deviation and survivorship fished
	  Nmta(fyr-lage+1,j) = Nfta(fyr-lage+1,j) * Imp(j); // use equilibrium male proportion at age to get initial number of males
    Nfta(fyr-lage+1,j) *= (1.-Imp(j)); // reduce number of females by males
	}
	

FUNCTION statedynamics
	
	// pop loop
	for(int i=fyr-lage+1;i<=lyr;i++)
	{
		Nfta(i+1,fage)=mfexp(log_rbar+wt(i));//Recruitment each year as rbar and deviation
		for(int a=fage;a<lage;a++)
		{
      Nfta(i+1,a+1)=Nfta(i,a)*S(i,a); // next year's females at age are are this year's surviving females at age
      Nmta(i+1,a+1)=Nmta(i,a)*S(i,a)+Nfta(i+1,a+1)*Pfm(a+1); // surviving males at age this year + females transitioning this year
      Nfta(i+1,a+1)*=(1.-Pfm(a+1));// adjust females
      rt_hat(i+1)=Nfta(i+1,fage); // observed recruitment
    }
		Nmta(i+1,lage)+=Nmta(i,lage)*S(i,lage); //Plus group calculation
	}
	
	// don't need summed pop metrics until fyr
	for(int i=fyr;i<=lyr;i++)
	{
    Nt(i)=Nfta(i)+Nmta(i);
    Nft(i)=sum(Nfta(i));
    Nmt(i)=sum(Nmta(i));
    ANft(i)=0.0;
    for(int a=fage_adult;a<=lage;a++){
	    ANft(i)+=Nfta(i,a);
	  }
  }
  Nta=Nfta+Nmta;
  ANt=ANft+Nmt;


FUNCTION fisheries_obs_model
	//cpue residuals ala waters and ludwig 1994
	rt_resid = log(Rt)-log(rt_hat(fyr-lage+1,lyr)); //Z values by year, assuming lognormal error around cpue
	q=mfexp(mean(rt_resid)); //q as exp(Z_bar), Ludwig and watlers z statistic
	rt_resid-=mean(rt_resid); //Z residuals as Z values by year - Zbar
	for(int i=fyr;i<=lyr;i++) pahat(i)=Nta(i)/sum(Nta(i));// proportions at age
	//cout<<pahat<<endl;
	//ad_exit(1);

FUNCTION ckmr_obs_model
	for(int i=fyr;i<=lyr;i++)
	{
	  Ef_denom(i)=0.0;
	  Em_denom(i)=0.0;
	  for(int a=fage;a<=lage;a++)
	  {
	    Rf(i,a)=Nfta(i,a)/Nft(i);
	    Rm(i,a)=Nmta(i,a)/Nmt(i);
	    Ef_denom(i)+=Rf(i,a)*Af(a)*Fec(a);
      Em_denom(i)+=Rm(i,a)*FecM(a);
	  }
	}
	

FUNCTION objective_function 
  dvar_vector p_vec(1,6);// penalties
  dvar_vector nll_vec(1,3);//vector for likelihood components
	nll_vec.initialize();
  
  // recruitment likelihood
  nll_vec(1)=dnorm(rt_resid,rsd)*index_weight; //residuals sd must be specified can not be estimated
  
  // comps likelihood
  double tau2;
	dvar_matrix nu(fyr,lyr,fage,lage);
	nll_vec(2)=dmvlogistic(comps,pahat,nu,tau2,0.001)*comps_weight;
	
	// ckmr likelihoods
	// cycle over the different cohort categories (i.e. the different values of c1 < c2)
  for(int i=fckmr;i<=nckmr;i++) {
    
    // for each cohort category, we have to cycle over the different possible values of a^*
    // to collect a sum.  We have some variables here for collecting those sums.
    Pffcc = 0.0;   // this is for P(\mathrm{HSP}_{f\rightarrow f}|c_1, c_2)
    Pfmcc = 0.0;
    Pmmcc = 0.0;
    int d = c2(i) - c1(i);
    
    // cycle over the different possible values of a^*.
    for(int as=d+1;as<=lage;as++) {

      // First deal with the mother-sharing half sibs
      if(Nft(c2(i)) * Rf(c2(i), as) > 0.0) {
        dvariable tmp_ff = 1.0;  // this is for accumulating a product that is P(\mathrm{HSP}_{f\rightarrow f}|c_1, c_2, a^*)
        tmp_ff *= Rf(c1(i), as-d) * Af(as-d) * Fec(as-d) / Ef_denom(c1(i));  // This is the E_{f, c_1, a^*-d} factor
        // then we have to cycle over the product of survival
        for(int s=0;s<d;s++) { tmp_ff *= S(c1(i)+s, as-d+s); }
        // then we have to cycle over the product of not transitioning to male between
        // giving birth initially and later.  We cycle from the time of the first-born's
        // birth PLUS one, to the time of the second-born's birth, because transitioning
        // happens before reproduction 
        for(int s=as-d+1;s<=as;s++) { tmp_ff *= (1.0 - Pfm(s)); }
        // then we multiply the 1 over the number of females at a^* term.
        tmp_ff *=  1.0 / (Nft(c2(i)) * Rf(c2(i), as));
        // and finally, add that to the running sum, while weighting by E_{f,c2,a^*}
        Pffcc += tmp_ff * Rf(c2(i), as) * Af(as) * Fec(as) / Ef_denom(c2(i));
      }
      
      // Now, deal with the half sibs in which the first-born's mother is the second-born's father
      if(Nmt(c2(i)) * Rm(c2(i), as) > 0.0) {
        dvariable tmp_fm = 1.0;  // this is for accumulating a product that is P(\mathrm{HSP}_{f\rightarrow f}|c_1, c_2, a^*)
        tmp_fm *= Rf(c1(i), as-d) * Af(as-d) * Fec(as-d) / Ef_denom(c1(i));  // This is the E_{f, c_1, a^*-d} factor
        // then we have to cycle over the product of survival
        for(int s=0;s<d;s++) { tmp_fm *= S(c1(i)+s, as-d+s); }
        // then we have to cycle over the product of not transitioning to male, but
        // multiply 1 minus that onto tmp_fm. We cycle from the time of the first-born's
        // birth PLUS one, to the time of the second-born's birth, because transitioning
        // happens before reproduction (when the individual is a year older than when it
        // produced the first-born)
        dvariable fmprod = 1.0;
        for(int s=as-d+1;s<=as;s++) { fmprod *= (1.0 - Pfm(s)); }
        tmp_fm *= (1.0 - fmprod);
        // then we multiply by the 1 over the number of males at a^* term.
        tmp_fm *=  1.0 / (Nmt(c2(i)) * Rm(c2(i), as));
        // and finally, add that to the running sum, while weighting by E_{m,c2,a^*}
        Pfmcc += tmp_fm * Rm(c2(i), as) * FecM(as) / Em_denom(c2(i));
      }
      
      // And now, deal with the half sibs in which the first-born's father is the second-born's father
      if(Nmt(c2(i)) * Rm(c2(i), as) > 0.0) {
        dvariable tmp_mm = 1.0;  // this is for accumulating a product that is P(\mathrm{HSP}_{f\rightarrow f}|c_1, c_2, a^*)
        tmp_mm *= Rm(c1(i), as-d) * FecM(as-d) / Em_denom(c1(i));  // This is the E_{m, c_1, a^*-d} factor
        // then we have to cycle over the product of survival
        for(int s=0;s<d;s++) { tmp_mm *= S(c1(i)+s, as-d+s); }
        // then we DO NOT have to cycle over anything about transitioning to male, because
        // once a male, always a male.
        // then we multiply by the 1 over the number of males at a^* term.
        tmp_mm *= 1.0 / (Nmt(c2(i)) * Rm(c2(i), as));
        // and finally, add that to the running sum, while weighting by E_{m,c2,a^*}
        Pmmcc += tmp_mm * Rm(c2(i), as) * FecM(as) / Em_denom(c2(i));
      }
    }
    
    // POPs
    Ppopffcc = 0.0;
    Ppopfmcc = 0.0;
    Ppopmmcc = 0.0;
    int dd = s1(i) - c2(i); // difference between sampling year of older inidividual and birth year of younger indiv
    
    if(dd = 0) // s1 sampled in c2's birth year
    {
      if(g1(i) == 0) // female at time of sampling
      {
        Ppopffcc = Fec(d) * Af(d) / (Ef_denom(c2(i)) * Nft(c2(i)) );
      } 
      else // male at time of sampling
      {
        Ppopmmcc = FecM(d) / (Em_denom(c2(i)) * Nmt(c2(i)) );
      }
    }
    else if(dd > 0) // s1 sampled after c2's birth year
    {
      if(g1(i) == 0) // female at time of sampling
      {
        // account for it not transitioning
        Ppopffcc = Fec(d) * Af(d) / (Ef_denom(c2(i)) * Nft(c2(i)) );
        // probability that it DIDN't transition from female to male between making offspring and being sampled
        for(int s=d+1;s<=s1(i)-c1(i);s++) { Ppopffcc *= (1.0 - Pfm(s)); }
      }
      else // male at time of sampling
      {
        // if male at time c2 -- don't incorporate prob of transitioning
        Ppopmmcc = FecM(d) / (Em_denom(c2(i)) * Nmt(c2(i)) );
        // if female at time c2 -- incorporate probability that it DID transition
        dvariable fmprod = 1.0;
        for(int s=d+1;s<=s1(i)-c1(i);s++) { fmprod *= (1.0 - Pfm(s)); }
        Ppopfmcc = (1.0 - fmprod) * Fec(d) * Af(d) / (Ef_denom(c2(i)) * Nft(c2(i)) );
      }
    }
    
    // Now, compute the HSD and HSS probs as per mtDNA heterozygosity
    PHSS = Pffcc + (1 - mtHZ_est) * (Pfmcc + Pmmcc);
    PHSD = mtHZ_est * (Pfmcc + Pmmcc);
    
    // compute PPHSS and PPHSD
    PPHSS = Ppopffcc + Ppopfmcc + (1 - mtHZ_est) * (Ppopmmcc);
    PPHSD = mtHZ_est * (Ppopmmcc);
    
    Pffcc_vec(i) = Pffcc;
    Pfmcc_vec(i) = Pfmcc;
    Pmmcc_vec(i) = Pmmcc;
    PHSS_vec(i) = PHSS;
    PHSD_vec(i) = PHSD;
    Ppopffcc_vec(i) = Ppopffcc;
    Ppopfmcc_vec(i) = Ppopfmcc;
    Ppopmmcc_vec(i) = Ppopmmcc;
    PPHSD_vec(i) = PPHSD;
    PPHSS_vec(i) = PPHSS;
    
    // and now we are ready to contribute to the negative log likelihood.  This is
    // just using the half siblings, currently.
    
    if(nll_type==0)
    {
      nll_vec(3) -=  log(PHSS + 1e-16) * HSS(i) + log(PHSD + 1e-16) * HSD(i) + log(1.0 - PHSS - PHSD) * (Not(i) + PO(i) + MO(i));
    } 
    else
    {
      nll_vec(3) -=  log(PPHSS + 1e-16) * MO(i) + log(PPHSD + 1e-16) * PO(i) + log(PHSS + 1e-16) * HSS(i) + log(PHSD + 1e-16) * HSD(i) + log(1.0 - PHSS - PHSD - PPHSS - PPHSD) * Not(i);
    }
  }
  nll_vec(3) *= ckmr_weight;
	
	p_vec.initialize();// create p_vec
	p_vec(1)=dlnorm(fbarc,-1.2,0.8); //
	p_vec(2)=dlnorm(fbarr,-1.2,0.8); //
	p_vec(3)=dlnorm(tsd,1.5,0.5); //
	
	// this keeps recruitment deviations constrained until in the last phase of estimation
	if(last_phase())
	{
		p_vec(4)=dnorm(wt,2);//week priors on the deviaitons
		p_vec(5)=dnorm(log_ftc_dev,2);
		p_vec(6)=dnorm(log_ftr_dev,2);
	}
	else
	{
		p_vec(4)=100.*norm2(wt); //stong priors preventing big deviations
		p_vec(5)=100.*norm2(log_ftc_dev);
		p_vec(6)=100.*norm2(log_ftr_dev);

	}
	
	nll=sum(nll_vec)+sum(p_vec);//total posterior as likelihoods and priors

FUNCTION mcmc_output
	if(iter==0)
	{
		ofstream ofs("MCMCoutput.mcmc");
		//ofs<<"Nft\t Nmt\t nll\t"<<endl;
	}
	iter++;
	ofstream ofs("MCMCoutput.mcmc",ios::app);
	ofs<<Nft<<"\t"<<Nmt<<"\t"<<nll<<endl;

REPORT_SECTION
	
  objective_function_value::pobjfun->gmax;
  REPORT(PHSS_vec);
  REPORT(PHSD_vec);
  REPORT(PPHSS_vec);
  REPORT(PPHSD_vec);
  REPORT(Pffcc_vec);
  REPORT(Pfmcc_vec);
  REPORT(Pmmcc_vec);
  REPORT(Ppopffcc_vec);
  REPORT(Ppopfmcc_vec);
  REPORT(Ppopmmcc_vec);
  REPORT(Ef_denom);
  REPORT(Nft);
  REPORT(Nmt);
  REPORT(Nfta);
  REPORT(Nmta);
  REPORT(lxf);
  REPORT(nll);
  REPORT(wt);
  REPORT(Rf);
  REPORT(Rm);
  REPORT(S);
  REPORT(Pfm);
  REPORT(Imp); 
  REPORT(FecM);
  REPORT(rt_hat);
  REPORT(q);
  REPORT(rt_resid);
  REPORT(Nta);
  REPORT(pahat);
  REPORT(comps);
  
  
TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <fvar.hpp>
	//#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;
	