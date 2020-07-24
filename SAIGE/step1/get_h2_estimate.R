load('casectrl.rda')
tau=modglmm$theta
tau
#[1] 1.0000000 0.6641714 #tau[1] is always 1 for binary traits

###for quantiative traits
h2 = tau[2]/(tau[1]+tau[2])


###binary traits
h2_liability=tau[2]/(tau[2]+pi^2/3)
h2_liability 
#[1] 0.1679729

###notes:
#the heritability is the point estimate for proportion of variance of the phenotype explained by the GRM
#h2 estimate for binary traits by SAIGE is underestimated 
#and the penalized quasi-likelihood used in SAIGE is known to be biased for heritability estimation but it works well for adjusting for sample-relatedness.
