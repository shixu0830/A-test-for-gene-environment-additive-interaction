################################################################
### Data application: ovarian cancer study
### Xu Shi (xushi@hsph.harvard.edu)
### Benedict H.W. Wong (wong01@g.harvard.edu)
### Tamar Sofer (tsofer@bwh.harvard.edu)
################################################################


rm(list=ls())
library(alr3)
library(maxLik)
library(MASS)
library(CGEN)
source("functions.R")
############### get data ###############
data(Xdata, package="CGEN")
dat = Xdata
dat$ethnicity.1=as.numeric(dat$ethnic.group == 1)
dat$ethnicity.2=as.numeric(dat$ethnic.group == 2)
dat$ethnicity.3=as.numeric(dat$ethnic.group == 3)
dat$famhist.1=as.numeric(dat$family.history == 1)
dat$famhist.2=as.numeric(dat$family.history == 2)
dat$age50=as.numeric(dat$age.group > 2) ## age at least 50
#### Adjusted covariates in the primary analysis
adjust.vars = c("age50","ethnicity.1","ethnicity.2","BRCA.history","famhist.1","famhist.2")
#### Adjusted covariates in the sensitivity analysis: remove BRCA.history
adjust.vars = c("age50","ethnicity.1","ethnicity.2","famhist.1","famhist.2")

rslt=NULL
for(children_or_oralcon in c(1,2)){
  for(ge.indep in c(T,F)){
    boot.dat=dat
    d = boot.dat$case.control
    g = boot.dat$BRCA.status
    if(children_or_oralcon==1){
      e = boot.dat$n.children
    }else{
      e = boot.dat$oral.years
    }
    X = cbind(intercept=1,as.matrix(boot.dat[,adjust.vars]))
    p=ncol(X)
    myg=g[d==0];mye=e[d==0];myX=as.matrix(X[d==0,]);
    reri.p=RERI.E.cont.X(d,g,e,X)
    temp <- additive.test(data = data.frame(d=d,g=g,e=as.numeric(e<mean(e)),boot.dat[,adjust.vars]), response.var = "d", snp.var="g", exposure.var="e", 
                          main.vars=adjust.vars,
                          op = list(indep = ge.indep))
    pvals.reri <- temp$RERI$pval
    pvals.han <- temp$pval.add
    
    
    ### our proposed tests
    if(ge.indep==T){
      ### assume independent
      MLE.rslt =  maxLik(loglik.indept,start=rep(0,p*2))
      par = MLE.rslt$estimate
      phi=0;alpha = par[1:( p )];beta = par[(1+p):( p*2 )]
    }else{
      ### not assume independent
      MLE.rslt = maxLik(loglik,start=rep(0,p*2+1))
      par = MLE.rslt$estimate
      phi=par[1];alpha = par[2:( 1+p )];beta = par[(2+p):( 1+p*2 )]
    }
    ### get u
    m.g = expit(X%*%alpha); m.e = exp(X%*%beta) 
    w=exp(-phi*g*e); u = w*(g - m.g)*(e - m.e)*d
    ###
    if(ge.indep==T){
      ### get IF
      d.u.d.theta = cbind(
        diag(as.numeric(  w*(e - m.e)*d*(-m.g*(1-m.g))  ))%*%X,
        diag(as.numeric(  w*(g - m.g)*d*(-m.e)  ))%*%X)
      E.d.u.d.theta = apply(d.u.d.theta,2,sum) ##order: (phi),g,e
      score=numDeriv::jacobian(loglik.indept,c(par))
    }else{
      ### get IF
      d.u.d.theta = cbind(-g*e*u,
                          diag(as.numeric(  w*(e - m.e)*d*(-m.g*(1-m.g))  ))%*%X,
                          diag(as.numeric(  w*(g - m.g)*d*(-m.e)  ))%*%X)
      E.d.u.d.theta = apply(d.u.d.theta,2,sum) ##order: phi,g,e
      score=numDeriv::jacobian(loglik,c(par))
    }
    score.IF=ginv(t(score))
    adj.IF = score.IF%*%E.d.u.d.theta
    u.IF = c(u[d==1],u[d==0]) + c(rep(0,length(which(d==1))),adj.IF)
    (p=(1 - pnorm( abs(  mean(u)*sqrt(length(u))/sd(u.IF)  )))*2)
    stat=sum(u)/sqrt(sum(u.IF^2)) 
    ( CI = round(mean(u)+c(-1,0,1)*qnorm(0.975)*sd(u.IF)/sqrt(length(u)),3) )
    
    
    rslt=rbind(rslt,c(sprintf("%.3f",phi),sprintf("%.3f",CI[2]),
                      paste0("(",sprintf("%.3f",CI[1]),",",sprintf("%.3f",CI[3]),")")
                      ,sprintf("%.3f",p),sprintf("%.3f",reri.p$pval)
                      ,sprintf("%.3f",pvals.reri),sprintf("%.3f",pvals.han)
    ))
  }
}
rslt
