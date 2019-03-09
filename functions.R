################################################################
### Functions
### Xu Shi (xushi@hsph.harvard.edu)
### Benedict H.W. Wong (wong01@g.harvard.edu)
### Tamar Sofer (tsofer@bwh.harvard.edu)
################################################################


expit = function(x){exp(x)/(1 + exp(x))}
logit = function(x){log(x/(1 - x))}


get.pval=function(g,e,X,d,ge.indep){
  if(ge.indep==T){
    ### assume G-E independence
    MLE.rslt =  maxLik(loglik.indept,start=rep(0,p*2))
    par = MLE.rslt$estimate
    phi=0;alpha = par[1:( p )];beta = par[(1+p):( p*2 )]
  }else{
    ### not assume G-E independence
    MLE.rslt = maxLik(loglik,start=rep(1,p*2+1))
    par = MLE.rslt$estimate
    phi=par[1];alpha = par[2:( 1+p )];beta = par[(2+p):( 1+p*2 )]
  }
  ### get u
  m.g = expit(X%*%alpha); m.e = exp(X%*%beta) 
  w=exp(-phi*g*e); u = w*(g - m.g)*(e - m.e)*d
  
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
  
  score.IF=score%*%solve(t(score)%*%(score)) #or score.IF=ginv(t(score))
  adj.IF = score.IF%*%E.d.u.d.theta
  u.IF = c(u[d==1],u[d==0]) + c(rep(0,length(which(d==1))),adj.IF)
  p=(1 - pnorm( abs(  mean(u)*sqrt(length(u))/sd(u.IF)  )))*2 #or p=(1 - pnorm( abs(  sum(u)/sqrt(sum(u.IF^2))  )))*2
  stat=sum(u)/sqrt(sum(u.IF^2)) 
  return(c(p=p,mean.u=mean(u),phi=phi,stat=stat,sd.u=sd(u.IF)/sqrt(length(u))
           ))
}


loglik.indept = function(par){#### Function for the log likelihood assuming G-E indep
  ##get coefs
  phi=0;alpha = par[1:( ncol(myX) )];beta = par[(1+ncol(myX)):( ncol(myX)*2 )]
  ##compute joint destribution f(G,E|X)
  lambdaX = exp(myX%*%beta)
  # ##proportional to f(G,E|X)
  fGE_cond_X = exp(myX%*%alpha*myg)/(1 + exp(myX%*%alpha))*  #fG_cond_E0X
    lambdaX^mye * exp(-lambdaX)/factorial(mye)*            #fE_cond_G0X
    exp(phi*myg*mye)                                       #OR(G,E|X), assume doesnt depend on X
  C_X = expit(myX%*%alpha)*exp( lambdaX*( exp(phi)-1 ) ) + 1/(1+ exp(myX%*%alpha)) ##normalizing const
  fGE_cond_X=fGE_cond_X/C_X
  return((log(fGE_cond_X)))
}
loglik = function(par){#### Function for the log likelihood without assuming G-E indep
  ##get coefs
  phi=par[1];alpha = par[2:( 1+ncol(myX) )];beta = par[(2+ncol(myX)):( 1+ncol(myX)*2 )]
  ##compute joint destribution f(G,E|X)
  lambdaX = exp(myX%*%beta)
  # ##proportional to f(G,E|X)
  fGE_cond_X = exp(myX%*%alpha*myg)/(1 + exp(myX%*%alpha))*#fG_cond_E0X
    lambdaX^mye * exp(-lambdaX)/factorial(mye)*            #fE_cond_G0X
    exp(phi*myg*mye)                                       #OR(G,E|X), assume doesnt depend on X
  C_X = expit(myX%*%alpha)*exp( lambdaX*( exp(phi)-1 ) ) + 1/(1+ exp(myX%*%alpha)) ##normalizing const
  fGE_cond_X=fGE_cond_X/C_X
  return((log(fGE_cond_X)))
}


RERI.E.cont = function(d,g,e){#### RERI using standard logistic regression for the outcome
  mod = glm(d ~ g + e + g*e, family = "binomial")
  temp = deltaMethod(mod, "exp(b1 + b2 + b3) - exp(b1) - exp(b2) + 1", parameterNames = paste("b", 0:3, sep = ""))
  test.stat = temp$Estimate/temp$SE
  pval = 2*(1 - pnorm(abs(test.stat)))
  return(list(test.stat = test.stat, pval = pval))
}
RERI.E.cont.X = function(d,g,e,X){#### RERI using standard logistic regression for the outcome with covariate adjustment
  X=X[,-1] #remove intercept
  colnames(X)=paste0("x",4:(ncol(X)+3))
  
  mod = glm(as.formula(paste("d~",paste(paste0("x",1:(ncol(X)+3)),collapse="+"))),data=data.frame(x1=g,x2=e,x3=g*e,X),family = "binomial")
  temp = deltaMethod(mod, "exp(x1+x2+x3) - exp(x1) - exp(x2) + 1", parameterNames = c(paste0("x",0:(ncol(X)+3))))
  test.stat = temp$Estimate/temp$SE
  pval = 2*(1 - pnorm(abs(test.stat)))
  return(list(test.stat = test.stat, pval = pval))
}

estimate.g.e.probs = function(g,e, ge.indep = F, return.equation = F, max.iter = 500, eps = 1e-6){
  #### Function for the proposed test under binary G and E (saturated model)
  converge = F
  iter = 0
  ind.alpha1 = 1
  ind.alpha2 = 2
  
  if (ge.indep){
    
    alpha.old = rep(0,2)
    
    U.vec = rep(0, length(alpha.old))
    U.deriv.mat = matrix(0, length(alpha.old), length(alpha.old))
    
    while (!converge & iter < max.iter){
      
      e1 = exp(alpha.old[ind.alpha1])
      e2 = exp(alpha.old[ind.alpha2])
      e12 = exp(sum(alpha.old))
      
      C = 1/(e12 + e1 + e2 + 1 ) 
      U.vec[ind.alpha1] = mean(g) - C*(e12 + e1)
      U.vec[ind.alpha2] = mean(e) - C*(e12 + e2)
      
      U.deriv.mat[ind.alpha1, ind.alpha1] = C^2*(e12 + e1)*(e12 + e1) - C*(e12 + e1)
      U.deriv.mat[ind.alpha1, ind.alpha2] = C^2*(e12 + e2)*(e12 + e1) - C*e12
      U.deriv.mat[ind.alpha2, ind.alpha1] = C^2*(e12 + e1)*(e12 + e2) - C*e12
      U.deriv.mat[ind.alpha2, ind.alpha2] = C^2*(e12 + e2)*(e12 + e2) - C*(e12 + e2)
      
      alpha.new = alpha.old - solve(U.deriv.mat) %*% U.vec
      
      if (max(abs(alpha.new - alpha.old)) <= eps) converge = T else{
        alpha.old = alpha.new
        iter = iter + 1
      }
      
    }
    
    
    if (return.equation){
      U = matrix(0, ncol = length(alpha.new), nrow = length(g))
      U[,ind.alpha1] = g - C*(exp(sum(alpha.new)) + exp(alpha.old[ind.alpha1]))
      U[,ind.alpha2] = e - C*(exp(sum(alpha.new)) + exp(alpha.old[ind.alpha2]))
      
      #		EU.inv.U = U %*% solve(U.deriv.mat) 
      EU.inv.U = U %*% solve(t(U) %*% U) 
      
    }
    
    alpha.new = c(alpha.new, 0)
    
  } else{
    
    alpha.old = rep(0,3)
    
    ind.alpha3 = 3
    
    U.vec = rep(0, length(alpha.old))
    U.deriv.mat = matrix(0, length(alpha.old), length(alpha.old))
    
    while (!converge & iter < max.iter){
      
      e1 = exp(alpha.old[ind.alpha1])
      e2 = exp(alpha.old[ind.alpha2])
      e123 = exp(sum(alpha.old))
      
      
      C = 1/(e123 + e1 + e2 + 1 ) 
      U.vec[ind.alpha1] = mean(g) - C*(e123 + e1)
      U.vec[ind.alpha2] = mean(e) - C*(e123 + e2)
      U.vec[ind.alpha3] = mean(g*e) - C*e123
      
      U.deriv.mat[ind.alpha1, ind.alpha1] = C^2*(e123 + e1)*(e123 + e1) - C*(e123 + e1)
      U.deriv.mat[ind.alpha1, ind.alpha2] = C^2*(e123 + e2)*(e123 + e1) - C*e123
      U.deriv.mat[ind.alpha1, ind.alpha3] = C^2*e123*(e123 + e1) - C*e123
      U.deriv.mat[ind.alpha2, ind.alpha1] = C^2*(e123 + e1)*(e123 + e2) - C*e123
      U.deriv.mat[ind.alpha2, ind.alpha2] = C^2*(e123 + e2)*(e123 + e2) - C*(e123 + e2)
      U.deriv.mat[ind.alpha2, ind.alpha3] = C^2*e123*(e123 + e2) - C*e123
      U.deriv.mat[ind.alpha3, ind.alpha1] = C^2*(e123 + e1)*e123 - C*e123
      U.deriv.mat[ind.alpha3, ind.alpha2] = C^2*(e123 + e2)*e123 - C*e123
      U.deriv.mat[ind.alpha3, ind.alpha3] = C^2*e123*e123 - C*e123
      
      alpha.new = alpha.old - solve(U.deriv.mat) %*% U.vec
      
      if (max(abs(alpha.new - alpha.old)) <= eps) converge = T else{
        alpha.old = alpha.new
        iter = iter + 1
      }
      
    }
    
    if (return.equation){
      U = matrix(0, ncol = length(alpha.new), nrow = length(g))
      U[,ind.alpha1] = g - C*(e123 + e1)
      U[,ind.alpha2] = e - C*(e123 + e2)
      U[,ind.alpha3] = g*e - C*e123		
      #		EU.inv.U = U %*% solve(U.deriv.mat) 
      EU.inv.U = U %*% solve(t(U) %*% U) 
      
    }
    
    
  }
  
  ##  joint probabilities:
  p.g1.e1 = C*exp(sum(alpha.new))
  p.g1.e0 = C*exp(alpha.new[ind.alpha1])
  p.g0.e1 = C*exp(alpha.new[ind.alpha2])
  p.g0.e0 = C
  
  ### return conditional probabilities:
  
  p.g1.cond.e1 = p.g1.e1/(p.g1.e1 + p.g0.e1)
  p.e1.cond.g1 = p.g1.e1/(p.g1.e0 + p.g1.e1)
  p.g1.cond.e0 = p.g1.e0/(p.g0.e0 + p.g1.e0)
  p.e1.cond.g0 = p.g0.e1/(p.g0.e1 + p.g0.e0)
  
  if (return.equation){
    return(list(cond.probs = list(p.g1.cond.e1 = p.g1.cond.e1, p.g1.cond.e0 = p.g1.cond.e0, p.e1.cond.g0 = p.e1.cond.g0, p.e1.cond.g1 = p.e1.cond.g1), coefs =  list(alpha = alpha.new) , equation = EU.inv.U))
  } else  {
    return(list(cond.probs = list(p.g1.cond.e1 = p.g1.cond.e1, p.g1.cond.e0 = p.g1.cond.e0, p.e1.cond.g0 = p.e1.cond.g0, p.e1.cond.g1 = p.e1.cond.g1), coefs = list(alpha = alpha.new) ))}
  
}