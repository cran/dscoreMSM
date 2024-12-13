#' Weibull baseline hazard
#' @param t time
#' @param shape shape parameter
#' @param scale scale parameter
#' @return hazard function value under Weibull distibution
weibulbh<-function(t,shape=2,scale=1) {
  shape<-shape
  scale<-scale
  return((shape/scale)*(t/scale)^(shape-1))
}

#' Exponential baseline hazard
#' @param t time
#' @param shape shape parameter
#' @return hazard function value under Exponential distibution
expbh<-function(t,shape=2){
  shape<-shape
  return(shape)
}

#' Gompartz baseline hazard
#' @param t time
#' @param shape shape parameter
#' @param scale scale parameter
#' @return hazard function value under Gompartz distibution
gompbh<-function(t,shape=2,scale=1){
  shape<-shape
  scale<-scale
  return(shape*exp(scale*t))
}

#' @title simulation of survival data
#' @description function for simulation of survival data assuming the data comes from a parametric coxph model with gamma frailty distribution
#' @param n number of individual
#' @param beta vector of regression coefficient for coxph model
#' @param fvar frailty variance value(currently the function works for gamma frailty only)
#' @param bhdist distribution of survival time at baseline e.g. "weibull","exponential","llogistic"
#' @param fdist distribution of frailty terms e.g. "gamma"
#' @param X model matrix for the coxPH model with particular choice of beta
#' @param ... user can assume the shape and scale parameter of baseline survival distribution
#' @details
#' The process for simulation of multistate survival data is described in our manuscript. As the process includes transition through different states and it involves simulating survival time in different transition. So we have demonstrated the code for simulation of simple survival model.
#' Suppose we want to simulate a survival data with parametric baseline hazard and parametric frailty model. The hazard model is as follows:
#'  \deqn{h_i(t)=z_ih_0(t)exp(\textbf{x}_i\beta)\;;i=1,2,3,...,n}
#' where the baseline survival time follow Weibull distribution and the hazard is
#' \deqn{h_0(t)=\rho \lambda t^{\rho-1}}. Similarly we can have Gompertz, log logistic distribution. The following are the formula for hazard and cummulative hazard function
#' For exponential: \eqn{h_0(t)=\lambda} and \eqn{H_0(t)=\lambda t}\;\eqn{\lambda>0}
#'     Gompertz: \eqn{h_0(t)=\lambda exp(\gamma t)} and \eqn{H_0(t)=\frac{\lambda}{\gamma}(exp(\gamma t)-1)};\eqn{\lambda,\gamma>0}
#'
#' @return simulated survival data for a single transition
#' @import stats
#' @export
#' @references
#' Vishwakarma, G. K., Bhattacherjee, A., Rajbongshi, B. K., & Tripathy, A. (2024). Censored imputation of time to event outcome through survival proximity score method. \emph{Journal of Computational and Applied Mathematics}, 116103;
#'
#' Bhattacharjee, A., Vishwakarma, G. K., Tripathy, A., & Rajbongshi, B. K. (2024). Competing risk multistate censored data modeling by propensity score matching method. \emph{Scientific Reports}, 14(1), 4368.
#' @examples
#' \donttest{
#' ##
#' n1<-1000
#' p1<-2
#' X1<-matrix(rnorm(n1*p1),n1,p1)
#' simulated_data<-simfdata(n=1000,beta=c(0.5,0.5),fvar=0.5,
#' X=X1)
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra K. Vishwakarma
#' @seealso \link{cphGM}
simfdata<- function(n,beta,fvar,bhdist='weibull',X,fdist='gamma',...) {
  p<-length(beta)
  lpred<-X%*%beta
  if(fdist!='gamma'){
    stop('currently this is the only distribution available for frailty')
  }
  if(!bhdist%in%c('weibull','exponential','llogistic')){
    stop('For baseline we have weibul,exponential,llogistic')
  }
  if(bhdist=='weibull'){
  frailty<-rgamma(n,shape=1/fvar,scale=fvar)
  cum_hazard<- function(t,i){
    exp(lpred[i])*frailty[i]*integrate(weibulbh,lower=0,upper=t)$value
  }
  s_times<-numeric(n)
  for(i in 1:n){
    U<-runif(1)
    s_times[i]<-uniroot(function(t)cum_hazard(t,i)+log(U),lower=0,upper=100)$root
  }
  c_times<-runif(n,min=0,max=max(s_times))
  observed_times<-pmin(s_times,c_times)
  status<-as.numeric(s_times<=c_times)
  data<-data.frame(id=1:n,X,frailty=frailty,time=observed_times,status=status)
  return(data)
  }else if(bhdist=='exponential'){
    frailty<-rgamma(n,shape=1/fvar,scale=fvar)
    cum_hazard<- function(t,i){
      exp(lpred[i])*frailty[i]*integrate(expbh,lower=0,upper=t)$value
    }
    s_times<-numeric(n)
    for(i in 1:n){
      U<-runif(1)
      s_times[i]<-uniroot(function(t)cum_hazard(t,i)+log(U),lower=0,upper=100)$root
    }
    c_times<-runif(n,min=0,max=max(s_times))
    observed_times<-pmin(s_times,c_times)
    status<-as.numeric(s_times<=c_times)
    data<-data.frame(id=1:n,X,frailty=frailty,time=observed_times,status=status)
    return(data)
  } else{
    frailty<-rgamma(n,shape=1/fvar,scale=fvar)
    cum_hazard<- function(t,i){
      exp(lpred[i])*frailty[i]*integrate(gompbh,lower=0,upper=t)$value
    }
    s_times<-numeric(n)
    for(i in 1:n){
      U<-runif(1)
      s_times[i]<-uniroot(function(t)cum_hazard(t,i)+log(U),lower=0,upper=100)$root
    }
    c_times<-runif(n,min=0,max=max(s_times))
    observed_times<-pmin(s_times,c_times)
    status<-as.numeric(s_times<=c_times)
    data<-data.frame(id=1:n,X,frailty=frailty,time=observed_times,status=status)
    return(data)
}
}
utils::globalVariables(c("rgamma","integrate","runif","uniroot"))

