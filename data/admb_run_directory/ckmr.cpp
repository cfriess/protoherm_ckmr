#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <ckmr.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
  datafile.allocate("datafile");
  ctlfile.allocate("ctlfile");
ad_comm::change_datafile_name(datafile);
  fyr.allocate("fyr");
  lyr.allocate("lyr");
  fage_adult.allocate("fage_adult");
  fage.allocate("fage");
  lage.allocate("lage");
  ages.allocate(fage,lage);
ages.fill_seqadd(fage,1); //create a sequence, fill vector starting with fage and incrementing by 1
  fckmr.allocate("fckmr");
  nckmr.allocate("nckmr");
  Strue.allocate(fyr-lage+1,lyr,fage,lage,"Strue");
  Af.allocate(fage,lage,"Af");
  Fec.allocate(fage,lage,"Fec");
  Rt.allocate(fyr-lage+1,lyr,"Rt");
  rsd.allocate("rsd");
  comps.allocate(fyr,lyr,fage,lage,"comps");
  vulr.allocate(fage,lage,"vulr");
  vulc.allocate(fage,lage,"vulc");
  Mage.allocate(fage,lage,"Mage");
  c1.allocate(fckmr,nckmr,"c1");
  c2.allocate(fckmr,nckmr,"c2");
  s1.allocate(fckmr,nckmr,"s1");
  g1.allocate(fckmr,nckmr,"g1");
  PO.allocate(fckmr,nckmr,"PO");
  MO.allocate(fckmr,nckmr,"MO");
  HSD.allocate(fckmr,nckmr,"HSD");
  HSS.allocate(fckmr,nckmr,"HSS");
  Not.allocate(fckmr,nckmr,"Not");
  mtHZ_est.allocate("mtHZ_est");
  et.allocate("et");
  ef.allocate("ef");
  comps_weight.allocate("comps_weight");
  index_weight.allocate("index_weight");
  ckmr_weight.allocate("ckmr_weight");
  nll_type.allocate("nll_type");
  eof.allocate("eof");
iter=0;
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
ad_comm::change_datafile_name(ctlfile);
  irbar.allocate("irbar");
  ifbarr.allocate("ifbarr");
  ifbarc.allocate("ifbarc");
  itinfl.allocate("itinfl");
  itsd.allocate("itsd");
  imfec_exp.allocate("imfec_exp");
  eofc.allocate("eofc");
		if(eofc!=999)
		{
			cout<<"Error reading control.\n Fix it."<<endl;
			ad_exit(1);
		}
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_rbar.allocate(6,17,1,"log_rbar");
  log_fbarc.allocate(-4.6,0.5,2,"log_fbarc");
  log_fbarr.allocate(-4.6,0.5,2,"log_fbarr");
  log_tinfl.allocate(1.5,3,-1,"log_tinfl");
  log_tsd.allocate(0.9,2.4,et,"log_tsd");
  log_mfec_exp.allocate(-4,2,ef,"log_mfec_exp");
  wt.allocate(fyr-lage,lyr,-10.,10.,3,"wt");
  log_ftc_dev.allocate(fyr-lage+1,lyr,-10.,10.,4,"log_ftc_dev");
  log_ftr_dev.allocate(fyr-lage+1,lyr,-10.,10.,4,"log_ftr_dev");
  nll.allocate("nll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
		log_rbar=log(irbar);
		log_fbarc=log(ifbarc);
		log_fbarr=log(ifbarr);
	  log_tinfl=log(itinfl);
		log_tsd=log(itsd);
		log_mfec_exp=log(imfec_exp);
  fbarc.allocate("fbarc");
  #ifndef NO_AD_INITIALIZE
  fbarc.initialize();
  #endif
  fbarr.allocate("fbarr");
  #ifndef NO_AD_INITIALIZE
  fbarr.initialize();
  #endif
  tinfl.allocate("tinfl");
  #ifndef NO_AD_INITIALIZE
  tinfl.initialize();
  #endif
  tsd.allocate("tsd");
  #ifndef NO_AD_INITIALIZE
  tsd.initialize();
  #endif
  mfec_exp.allocate("mfec_exp");
  #ifndef NO_AD_INITIALIZE
  mfec_exp.initialize();
  #endif
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
  #endif
  lxf.allocate(fage,lage,"lxf");
  #ifndef NO_AD_INITIALIZE
    lxf.initialize();
  #endif
  ftc.allocate(fyr-lage+1,lyr,"ftc");
  #ifndef NO_AD_INITIALIZE
    ftc.initialize();
  #endif
  ftr.allocate(fyr-lage+1,lyr,"ftr");
  #ifndef NO_AD_INITIALIZE
    ftr.initialize();
  #endif
  ft.allocate(fyr-lage+1,lyr,"ft");
  #ifndef NO_AD_INITIALIZE
    ft.initialize();
  #endif
  Z.allocate(fyr-lage+1,lyr,fage,lage,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(fyr-lage+1,lyr,fage,lage,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  Frec.allocate(fyr-lage+1,lyr,fage,lage,"Frec");
  #ifndef NO_AD_INITIALIZE
    Frec.initialize();
  #endif
  Fcom.allocate(fyr-lage+1,lyr,fage,lage,"Fcom");
  #ifndef NO_AD_INITIALIZE
    Fcom.initialize();
  #endif
  Nfta.allocate(fyr-lage+1,lyr+1,fage,lage,"Nfta");
  #ifndef NO_AD_INITIALIZE
    Nfta.initialize();
  #endif
  Nmta.allocate(fyr-lage+1,lyr+1,fage,lage,"Nmta");
  #ifndef NO_AD_INITIALIZE
    Nmta.initialize();
  #endif
  Nta.allocate(fyr-lage+1,lyr+1,fage,lage,"Nta");
  #ifndef NO_AD_INITIALIZE
    Nta.initialize();
  #endif
  Nt.allocate(fyr,lyr+1,fage,lage,"Nt");
  #ifndef NO_AD_INITIALIZE
    Nt.initialize();
  #endif
  Nft.allocate(fyr,lyr+1,"Nft");
  Nmt.allocate(fyr,lyr+1,"Nmt");
  ANft.allocate(fyr,lyr+1,"ANft");
  ANt.allocate(fyr,lyr+1,"ANt");
  Rf.allocate(fyr,lyr,fage,lage,"Rf");
  #ifndef NO_AD_INITIALIZE
    Rf.initialize();
  #endif
  Rm.allocate(fyr,lyr,fage,lage,"Rm");
  #ifndef NO_AD_INITIALIZE
    Rm.initialize();
  #endif
  Pfm.allocate(fage,lage,"Pfm");
  #ifndef NO_AD_INITIALIZE
    Pfm.initialize();
  #endif
  Imp.allocate(fage,lage,"Imp");
  #ifndef NO_AD_INITIALIZE
    Imp.initialize();
  #endif
  FecM.allocate(fage,lage,"FecM");
  #ifndef NO_AD_INITIALIZE
    FecM.initialize();
  #endif
  NF.allocate(fage,lage,"NF");
  #ifndef NO_AD_INITIALIZE
    NF.initialize();
  #endif
  NM.allocate(fage,lage,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  pahat.allocate(fyr,lyr,fage,lage,"pahat");
  #ifndef NO_AD_INITIALIZE
    pahat.initialize();
  #endif
  Pffcc.allocate("Pffcc");
  #ifndef NO_AD_INITIALIZE
  Pffcc.initialize();
  #endif
  Pfmcc.allocate("Pfmcc");
  #ifndef NO_AD_INITIALIZE
  Pfmcc.initialize();
  #endif
  Pmmcc.allocate("Pmmcc");
  #ifndef NO_AD_INITIALIZE
  Pmmcc.initialize();
  #endif
  Ppopffcc.allocate("Ppopffcc");
  #ifndef NO_AD_INITIALIZE
  Ppopffcc.initialize();
  #endif
  Ppopfmcc.allocate("Ppopfmcc");
  #ifndef NO_AD_INITIALIZE
  Ppopfmcc.initialize();
  #endif
  Ppopmmcc.allocate("Ppopmmcc");
  #ifndef NO_AD_INITIALIZE
  Ppopmmcc.initialize();
  #endif
  PHSS.allocate("PHSS");
  #ifndef NO_AD_INITIALIZE
  PHSS.initialize();
  #endif
  PHSD.allocate("PHSD");
  #ifndef NO_AD_INITIALIZE
  PHSD.initialize();
  #endif
  PPHSS.allocate("PPHSS");
  #ifndef NO_AD_INITIALIZE
  PPHSS.initialize();
  #endif
  PPHSD.allocate("PPHSD");
  #ifndef NO_AD_INITIALIZE
  PPHSD.initialize();
  #endif
  Ef_denom.allocate(fyr,lyr,"Ef_denom");
  #ifndef NO_AD_INITIALIZE
    Ef_denom.initialize();
  #endif
  Em_denom.allocate(fyr,lyr,"Em_denom");
  #ifndef NO_AD_INITIALIZE
    Em_denom.initialize();
  #endif
  Pffcc_vec.allocate(fckmr,nckmr,"Pffcc_vec");
  #ifndef NO_AD_INITIALIZE
    Pffcc_vec.initialize();
  #endif
  Pfmcc_vec.allocate(fckmr,nckmr,"Pfmcc_vec");
  #ifndef NO_AD_INITIALIZE
    Pfmcc_vec.initialize();
  #endif
  Pmmcc_vec.allocate(fckmr,nckmr,"Pmmcc_vec");
  #ifndef NO_AD_INITIALIZE
    Pmmcc_vec.initialize();
  #endif
  PHSS_vec.allocate(fckmr,nckmr,"PHSS_vec");
  #ifndef NO_AD_INITIALIZE
    PHSS_vec.initialize();
  #endif
  PHSD_vec.allocate(fckmr,nckmr,"PHSD_vec");
  #ifndef NO_AD_INITIALIZE
    PHSD_vec.initialize();
  #endif
  Ppopffcc_vec.allocate(fckmr,nckmr,"Ppopffcc_vec");
  #ifndef NO_AD_INITIALIZE
    Ppopffcc_vec.initialize();
  #endif
  Ppopfmcc_vec.allocate(fckmr,nckmr,"Ppopfmcc_vec");
  #ifndef NO_AD_INITIALIZE
    Ppopfmcc_vec.initialize();
  #endif
  Ppopmmcc_vec.allocate(fckmr,nckmr,"Ppopmmcc_vec");
  #ifndef NO_AD_INITIALIZE
    Ppopmmcc_vec.initialize();
  #endif
  PPHSD_vec.allocate(fckmr,nckmr,"PPHSD_vec");
  #ifndef NO_AD_INITIALIZE
    PPHSD_vec.initialize();
  #endif
  PPHSS_vec.allocate(fckmr,nckmr,"PPHSS_vec");
  #ifndef NO_AD_INITIALIZE
    PPHSS_vec.initialize();
  #endif
  rt_hat.allocate(fyr-lage+1,lyr+1,"rt_hat");
  rt_resid.allocate(fyr-lage+1,lyr,"rt_resid");
  #ifndef NO_AD_INITIALIZE
    rt_resid.initialize();
  #endif
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

void model_parameters::userfunction(void)
{
  nll =0.0;
	initialization();
	statedynamics();
	fisheries_obs_model();
	ckmr_obs_model();
	objective_function();
	if(mceval_phase())
	{ 
		mcmc_output();
	}
}

void model_parameters::initialization(void)
{
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
}

void model_parameters::statedynamics(void)
{
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
}

void model_parameters::fisheries_obs_model(void)
{
	//cpue residuals ala waters and ludwig 1994
	rt_resid = log(Rt)-log(rt_hat(fyr-lage+1,lyr)); //Z values by year, assuming lognormal error around cpue
	q=mfexp(mean(rt_resid)); //q as exp(Z_bar), Ludwig and watlers z statistic
	rt_resid-=mean(rt_resid); //Z residuals as Z values by year - Zbar
	for(int i=fyr;i<=lyr;i++) pahat(i)=Nta(i)/sum(Nta(i));// proportions at age
	//cout<<pahat<<endl;
	//ad_exit(1);
}

void model_parameters::ckmr_obs_model(void)
{
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
}

void model_parameters::objective_function(void)
{
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
}

void model_parameters::mcmc_output(void)
{
	if(iter==0)
	{
		ofstream ofs("MCMCoutput.mcmc");
		//ofs<<"Nft\t Nmt\t nll\t"<<endl;
	}
	iter++;
	ofstream ofs("MCMCoutput.mcmc",ios::app);
	ofs<<Nft<<"\t"<<Nmt<<"\t"<<nll<<endl;
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
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
}

void model_parameters::final_calcs()
{
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
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
