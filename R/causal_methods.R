weighted_mean <- function(x,w) sum(x*w)/sum(w)
smd <- function(x,a,w=rep(1,length(a))) {
  m1<-weighted_mean(x[a==1],w[a==1]);m0<-weighted_mean(x[a==0],w[a==0])
  v1<-weighted_mean((x[a==1]-m1)^2,w[a==1]);v0<-weighted_mean((x[a==0]-m0)^2,w[a==0])
  (m1-m0)/sqrt((v1+v0)/2)
}

fit_propensity <- function(d) {
  fit<-glm(treatment~age+female+comorbidity+severity+prior_event+healthcare_use,
           family=binomial(),data=d)
  pmin(pmax(predict(fit,type="response"),.01),.99)
}

estimate_effects <- function(d) {
  n<-nrow(d); a<-d$treatment; y<-d$outcome; ps<-fit_propensity(d)
  p_a<-mean(a); sw<-ifelse(a==1,p_a/ps,(1-p_a)/(1-ps)); cap<-quantile(sw,c(.01,.99));sw<-pmin(pmax(sw,cap[1]),cap[2])
  crude_if<-a*(y-mean(y[a==1]))/mean(a)-(1-a)*(y-mean(y[a==0]))/mean(1-a)
  crude<-mean(y[a==1])-mean(y[a==0])
  outcome_fit<-glm(outcome~treatment+age+female+comorbidity+severity+prior_event+healthcare_use,family=binomial(),data=d)
  d1<-d0<-d;d1$treatment<-1;d0$treatment<-0;m1<-predict(outcome_fit,d1,type="response");m0<-predict(outcome_fit,d0,type="response")
  reg<-mean(m1-m0)
  x1<-model.matrix(delete.response(terms(outcome_fit)),d1);x0<-model.matrix(delete.response(terms(outcome_fit)),d0)
  gradient<-colMeans(m1*(1-m1)*x1-m0*(1-m0)*x0)
  reg_se<-sqrt(drop(t(gradient)%*%vcov(outcome_fit)%*%gradient))
  iptw<-weighted_mean(y[a==1],sw[a==1])-weighted_mean(y[a==0],sw[a==0])
  iptw_if<-a*sw*(y-weighted_mean(y[a==1],sw[a==1]))/mean(a*sw)-(1-a)*sw*(y-weighted_mean(y[a==0],sw[a==0]))/mean((1-a)*sw)
  dr_score<-m1-m0+a*(y-m1)/ps-(1-a)*(y-m0)/(1-ps);dr<-mean(dr_score)
  controls<-which(a==0);treated<-which(a==1);ord<-controls[order(ps[controls])];pos<-findInterval(ps[treated],ps[ord]);pos<-pmax(1,pmin(length(ord),pos));pos2<-pmin(length(ord),pos+1);choose2<-abs(ps[ord[pos2]]-ps[treated])<abs(ps[ord[pos]]-ps[treated]);match_control<-ord[ifelse(choose2,pos2,pos)];diff<-y[treated]-y[match_control]
  vals<-data.frame(method=c("Crude","Regression adjusted","Propensity score matching (ATT)","Stabilized IPTW (ATE)","Doubly robust AIPW (ATE)"),estimate=c(crude,reg,mean(diff),iptw,dr),se=c(sd(crude_if)/sqrt(n),reg_se,sd(diff)/sqrt(length(diff)),sd(iptw_if)/sqrt(n),sd(dr_score)/sqrt(n)))
  vals$lower<-vals$estimate-1.96*vals$se;vals$upper<-vals$estimate+1.96*vals$se
  list(effects=vals,ps=ps,weights=sw,true_ate=attr(d,"true_ate"))
}

balance_table <- function(d, ps, weights) {
  vars<-c("age","female","comorbidity","severity","prior_event","healthcare_use")
  data.frame(variable=vars,before=sapply(d[vars],smd,a=d$treatment),after=sapply(d[vars],smd,a=d$treatment,w=weights),row.names=NULL)
}
