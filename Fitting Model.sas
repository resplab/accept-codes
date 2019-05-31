%LET  B_LIN_F_TEMPLATE = b0
	                    + b_male*gender
	                    + b_age10*age10
	                    + b_nowsmk*nowsmk
	                    + b_oxygen*oxygen
	                    + b_fev1pp100*fev1pp100
	                    + b_sgrq10*sgrq10
	                    + b_BMI10*BMI10
	                    + b_statin*indicated_statin
	                    + b_randomized_azithromycin*randomized_azithromycin
	                    + b_LAMA*LAMA
	                    + b_LABA*LABA
	                    + b_ICS*ICS
						+ b_randomized_LAMA*randomized_LAMA
	                    + b_randomized_LABA*randomized_LABA
	                    + b_randomized_ICS*randomized_ICS
						+ b_randomized_statin*randomized_statin;

%LET C_LIN_F_TEMPLATE = c0
	                    + c_male*gender
	                    + c_age10*age10
	                    + c_nowsmk*nowsmk
	                    + c_oxygen*oxygen
	                    + c_fev1pp100*fev1pp100
	                    + c_BMI10*BMI10
	                    + c_sgrq10*sgrq10
	                    + c_statin*indicated_statin
	                    + c_randomized_azithromycin*randomized_azithromycin
	                    + c_LAMA*LAMA
	                    + c_LABA*LABA
	                    + c_ICS*ICS
					    + c_randomized_LAMA*randomized_LAMA
	                    + c_randomized_LABA*randomized_LABA
	                    + c_randomized_ICS*randomized_ICS
						+ c_randomized_statin*randomized_statin;

TITLE "Full Random Effect Model" ;
PROC NLMIXED DATA=exacEvents empirical cov;
   *initial values from fixed effect;
   PARMS gamma=0.925
            b0=0.2573 b_male=-0.06741 b_age10=-0.01959 b_nowsmk=-0.2196 b_oxygen=0.0165 b_fev1pp100=-0.41 b_sgrq10=0.1022
     		b_statin=0.1382 b_randomized_azithromycin=-0.1047 b_LAMA=0.1011 b_LABA=0.08752 b_ICS=0.1881 b_BMI10=-0.1043
			b_randomized_LAMA=0.1741 b_randomized_LABA=0.1175 b_randomized_ICS=-0.2792 b_randomized_statin=-0.09399 
            
            c0=-1.9209 c_male=0.4362 c_age10=-0.01768 c_nowsmk=0.1239 c_oxygen=0.3602 c_fev1pp100=-0.3644 c_sgrq10=0.15
            c_statin=0.227 c_randomized_azithromycin=-0.1534 c_LAMA=-0.1682 c_LABA=-0.1932 c_ICS=0.2154 c_BMI10=-0.04688
     		c_randomized_LAMA=0.06467 c_randomized_LABA=-0.3852 c_randomized_ICS=-0.01927 c_randomized_statin=0.06780  

            v1=1 v2=1 cov=0;  

	BOUNDS gamma>0,v1>0, v2>0;

    ***RATE component***;
    *b_lin_f: the fixed-effects linear predictor for rate component;
     b_lin_f = &B_LIN_F_TEMPLATE;

    *z1 is the random-effect term for rate;
    b_lin=b_lin_f + z1;
    alpha = exp(b_lin);

    *Survival and hazard functions of Weibull;
	S_t = EXP(-(alpha * tte0) ** gamma);
    h_t = gamma*alpha*((tte0*alpha)**(gamma-1));

    *ll1: the log-likelihood for the rate component;
    ll1 = (event>0)*log(h_t)+(event=0)*log(S_t);

    ***SEVERITY component***;
    *c_lin_f: the fixed-effects linear predictor for severity component;
    c_lin_f = &C_LIN_F_TEMPLATE; 

    *z2 is the random-effect term for severity;
    c_lin=c_lin_f + z2;

    *p is the probability of having severe exacerbation, conditional on having an exacerbation;
    p=EXP(c_lin)/(1+EXP(c_lin));
        *ll2: the log-likelihood for the severity component;
        *If for the censoring even (coded as event=0), there is no contribution to the severity component;
    ll2=LOG((event=0)*1+(event>0)*((event<=1)*(1-p)+(event>1)*p));

    ***Likelihood calculations***;
    ll=ll1+ll2;
    *We use "general" as we have constructed the likelihood outselves;
    model event ~ general(ll);

	random z1 z2 ~ normal([0,0],[v1,cov,v2]) subject=id OUT=_z;

    ****Log Hazard Ratios for Rate Component ****;
    ESTIMATE "ln(HR_male)" (b_male*gamma);
    ESTIMATE "ln(HR_age)" (b_age10*gamma);
    ESTIMATE "ln(HR_nowmsk)" (b_nowsmk*gamma);
    ESTIMATE "ln(HR_oxygen)" (b_oxygen*gamma);
    ESTIMATE "ln(HR_fev1pp100)" (b_fev1pp100*gamma);
    ESTIMATE "ln(HR_sgrq10)" (b_sgrq10*gamma);
    ESTIMATE "ln(HR_BMI10)" (b_BMI10*gamma);
    ESTIMATE "ln(HR_statin)" (b_statin*gamma);
    ESTIMATE "ln(HR_LAMA)" (b_LAMA*gamma);
    ESTIMATE "ln(HR_LABA)" (b_LABA*gamma);
    ESTIMATE "ln(HR_ICS)" (b_ICS*gamma);

    ****Log Odds Ratios for Severity Component ****;
    ESTIMATE "ln(OR_male)" (c_male*gamma);
    ESTIMATE "ln(OR_age)" (c_age10*gamma);
    ESTIMATE "ln(OR_nowmsk)" (c_nowsmk*gamma);
    ESTIMATE "ln(OR_oxygen)" (c_oxygen*gamma);
    ESTIMATE "ln(OR_fev1pp100)" (c_fev1pp100*gamma);
    ESTIMATE "ln(OR_sgrq10)" (c_sgrq10*gamma);
    ESTIMATE "ln(OR_BMI10)" (c_BMI10*gamma);
    ESTIMATE "ln(OR_statin)" (c_statin*gamma);
    ESTIMATE "ln(OR_LAMA)" (c_LAMA*gamma);
    ESTIMATE "ln(OR_LABA)" (c_LABA*gamma);
    ESTIMATE "ln(OR_ICS)" (c_ICS*gamma);

RUN;
title; 
