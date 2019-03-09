################################################################
### Simulation code for Binary G and Count E:
### Xu Shi (xushi@hsph.harvard.edu)
### Benedict H.W. Wong (wong01@g.harvard.edu)
### Tamar Sofer (tsofer@bwh.harvard.edu)
################################################################

rm(list=ls())
library(alr3)
library(maxLik)
library(MASS)
library(boot)
require(CGEN)
source("functions.R")
a0.seq = c(logit(0.01),logit(0.05),logit(0.1))  ## baseline risk
g.mean.seq=c(0.5,0.2,0.05)
signal.seq = seq(0, 0.5, by = 0.1) ##0:typeIerror; >0:power 
## we here look at a few combinations of gene/environement effects
env.eff.seq = c(log(0.7), log(1.2), log(2))  
gene.eff.seq = c(log(0.7), log(1.2), log(2))
n.sim=10000
nn=4000
N = nn/0.001
n.cont = nn
n.case = nn
n=n.cont*2
pval=pval2=pval.indept=NULL
for(i in 1:n.sim){
  #### simulate data
  set.seed((myseed-1)*n.sim+i)
  e = rpois(N,lambda=1)
  g = rbinom(N,1,g.mean)
  pD = expit(a0+gene.eff*g)+expit(a0+env.eff*e)+expit(a0)
  D = rbinom(length(pD),1,pD)
  cases = as.data.frame(cbind(e,g,D)[sample(which(D == 1), n.case),])
  conts = as.data.frame(cbind(e,g,D)[sample(which(D == 0), n.cont),])
  d = c(cases$D, conts$D)
  e = c(cases$e, conts$e)
  g = c(cases$g, conts$g)
  X = matrix(rep(1,length(g)),nrow=length(g))
  p=ncol(X)
  myg=g[d==0];mye=e[d==0];myX=as.matrix(X[d==0,]);
  
  #### Proposed test without assuming G-E independence
  pval.i = get.pval(g,e,X,d,ge.indep=F)
  pval=c(pval,pval.i["p"])
  #### Proposed test assuming G-E independence
  pval.indept.i = get.pval(g,e,X,d,ge.indep=T)
  if(inherits(pval.indept.i, "try-error")){pval.indept.i=NA}
  pval.indept=c(pval.indept,pval.indept.i["p"])
  #### RERI using standard logistic regression
  reri=RERI.E.cont(d,g,e)
  pval2=c(pval2,reri$pval)
}

#### compute type I error at sig level of 0.05
mean(pval<0.05)
mean(pval.indept<0.05)
mean(pval2<0.05)

