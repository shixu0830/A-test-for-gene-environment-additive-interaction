################################################################
### Simulation code for Binary G and Binary E:
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
### select one comb of the above seq
gene.eff = gene.eff.seq[1]  #1,2,3
env.eff = env.eff.seq[1]    #1,2,3
a0 = a0.seq[1]              #1,2
g.mean = g.mean.seq[1]      #1,2,3
sig = signal.seq[1]         #1,2,3
## setting the coeffcient of g*e such that RERI = sig.
gene.env.eff = logit( (sig -1)*expit(a0) + expit(a0 + env.eff) + expit(a0 + gene.eff)) - a0 - gene.eff - env.eff  
n.sim=10000
nn=4000
N <- nn/0.001
n.cont <- nn
n.case <- nn
n=n.cont*2

test.stat.U = pval.U=test.stat.U.ind = pval.U.ind=
  pvals.reri=pvals.UML=pvals.mult=pvals.wald.mult.indep=pvals.mult.indep=
  pvals.han=pvals.han.indep=test.stat.RERI=pval.RERI=
  rep(NA, n.sim)

for(i in 1:n.sim){
  set.seed(i)
  e <- rbinom(N,1,0.2)
  g <- rbinom(N,1,g.mean)
  pD <- expit(a0+env.eff*e+gene.eff*g+gene.env.eff*g*e)
  D <- rbinom(length(pD),1,pD)
  cases <- as.data.frame(cbind(e,g,D)[sample(which(D == 1), n.case),])
  conts <- as.data.frame(cbind(e,g,D)[sample(which(D == 0), n.cont),])
  d <- c(cases$D, conts$D)
  e <- c(cases$e, conts$e)
  g <- c(cases$g, conts$g)
  
  ############################################################
  ### Han and RERI
  temp <- additive.test(data = data.frame(d=d,g=g,e=e), response.var = "d", snp.var="g", exposure.var="e")
  temp.indep <- additive.test(data = data.frame(d=d,g=g,e=e), response.var = "d", snp.var="g", exposure.var="e", op = list(indep = T))
  pvals.reri[i] <- temp$RERI$pval
  pvals.han[i] <- temp$pval.add 
  pvals.han.indep[i] <- temp.indep$pval.add
  
  ############################################################
  ## RERI
  temp.RERI <- RERI.E.cont (d,g,e)
  test.stat.RERI[i] <- temp.RERI$test.stat
  pval.RERI[i] <- temp.RERI$pval
  
  ############################################################
  ## Proposed test stat U
  ### first: test that assumes GE independence
  probs.mod <- estimate.g.e.probs(g[which(d == 0)], e[which(d == 0)], ge.indep = T, return.equation = T)
  p10 <- probs.mod$cond.probs$p.g1.cond.e0
  p11 <- probs.mod$cond.probs$p.g1.cond.e1
  p20 <- probs.mod$cond.probs$p.e1.cond.g0
  
  u <- (g-p10)*(e-p20)*d
  mean.u = mean(u)
  alpha <- probs.mod$coef$alpha   ## alpha is the parameter of the model f(g,e|D=0) 
  exp.2 <- exp(alpha[2])
  exp.1 <- exp(alpha[1])
  ## derivative of U w.r.t. alpha
  partial.u.partial.V.p <- cbind(
    d*(0 - exp.1/((exp.1 + 1)^2) )*(e - exp.2/(1+exp.2)),
    d*(g - exp.1/(exp.1 + 1))*(0 - exp.2/(( 1+ exp.2)^2))
  )
  # its mean
  E.partial.u.partial.V.p <-  colMeans(partial.u.partial.V.p)
  alpha <- alpha[1:2]
  ## the equations accounting for estimation of alpha
  term.account.p <-  n*rbind(matrix(0, nrow = sum(d == 1), ncol = length(alpha) ), probs.mod$equation) %*% E.partial.u.partial.V.p	
  u.account.p.ind <- u - term.account.p
  sd.u <- sqrt(mean(u.account.p.ind^2)/(n))
  test.stat.U.ind[i] <- mean.u/sd.u
  pval.U.ind[i]=(1 - pnorm( abs(mean.u/sd.u)))*2
  ### second: test that does not assumes GE independence
  probs.mod <- estimate.g.e.probs(g[which(d == 0)], e[which(d == 0)], ge.indep = F, return.equation = T)
  p10 <- probs.mod$cond.probs$p.g1.cond.e0
  p11 <- probs.mod$cond.probs$p.g1.cond.e1
  p20 <- probs.mod$cond.probs$p.e1.cond.g0
  alpha <- probs.mod$coef$alpha
  u <- exp(-alpha[3]*g*e)*(g-p10)*(e-p20)*d
  mean.u = mean(u)
  exp.2 <- exp(alpha[2])
  exp.1 <- exp(alpha[1])
  term.3 <- exp(-alpha[3]*g*e)
  ## derivative of U w.r.t. alpha
  partial.u.partial.V.p <- cbind(
    d*term.3*(0 - exp.1/((exp.1 + 1)^2) )*(e - exp.2/(1+exp.2)),
    d*term.3*(g - exp.1/(exp.1 + 1))*(0 - exp.2/(( 1+ exp.2)^2)),
    -d*term.3*g*e*(g - exp.1/(exp.1 + 1))*(e - exp.2/(1+exp.2))
  )
  # its mean
  E.partial.u.partial.V.p <-  colMeans(partial.u.partial.V.p)
  term.account.p <-  n*rbind(matrix(0, nrow = sum(d == 1), ncol = length(alpha) ), probs.mod$equation) %*% E.partial.u.partial.V.p	
  u.account.p <- u - term.account.p
  sd.u <- sqrt(mean(u.account.p^2)/n)
  test.stat.U[i] <- mean.u/sd.u
  pval.U[i]=(1 - pnorm( abs(mean.u/sd.u)))*2
  if(i%%(n.sim/10)==0){print(i)}
}


