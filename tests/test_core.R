source("R/simulate_cohort.R");source("R/causal_methods.R")
d<-simulate_rwe_cohort(n=2000);a<-estimate_effects(d);b<-balance_table(d,a$ps,a$weights)
stopifnot(nrow(d)==2000,nrow(a$effects)==5,all(a$ps>0&a$ps<1),all(is.finite(a$effects$estimate)),max(abs(b$after))<max(abs(b$before)),attr(d,"true_ate")<0)
