#' @title CoxPH model with parametric baseline and frailty terms
#' @description Function for estimating the parameters of coxPH model with frailty terms
#' @param formula survival model formula like Surv(time,status)~x1+x2
#' @param fterm frailty term like c('gamma','center'). Currently we have the option for gamma distribution.
#' @param bhdist distribution of survival time at baseline. Available option 'weibull','exponential','gompertz',
#' @param Time survival time column
#' @param status survival status column
#' @param id id column
#' @param data dataset
#' @param method options are 'LFGS','L-BFGS-G','CG' etc. for more details see \link{optim}
#' @param maxit maximum number of iteration
#' @details
#' The hazard model is as follows:
#'  \deqn{h_i(t)=z_ih_0(t)exp(\textbf{x}_i\beta)\;;i=1,2,3,...,n}
#' where baseline survival distribution could be Weibull distribution and the hazard function is:
#' \deqn{h_0(t)=\rho \lambda t^{\rho-1}}. Similarly we can have Expoenetial, log logistic distribution. The following are the formula for hazard and cummulative hazard function
#' For exponential: \eqn{h_0(t)=\lambda} and \eqn{H_0(t)=\lambda t}\;\eqn{\lambda>0}
#'     Gompertz: \eqn{h_0(t)=\lambda exp(\gamma t)} and \eqn{H_0(t)=\frac{\lambda}{\gamma}(exp(\gamma t)-1)};\eqn{\lambda,\gamma>0}
#' The frailty term \eqn{z_i} follows Gamma distribution with parameter \eqn{\theta}. The parameter estimates are obtained by maximising the log likelihood
#' \deqn{\prod_{i=1}^{n}l_i(\beta,\theta,\lambda,\rho)}
#' The method argument allows the user to select suitable optimisation method available in \code{optim} function.
#' @return Estimates obtained from coxph model with the frailty terms.
#' @import stats
#' @export
#' @references
#' Vishwakarma, G. K., Bhattacherjee, A., Rajbongshi, B. K., & Tripathy, A. (2024). Censored imputation of time to event outcome through survival proximity score method. \emph{Journal of Computational and Applied Mathematics}, 116103;
#'
#' Bhattacharjee, A., Vishwakarma, G. K., Tripathy, A., & Rajbongshi, B. K. (2024). Competing risk multistate censored data modeling by propensity score matching method. \emph{Scientific Reports}, 14(1), 4368.
#' @examples
#' \donttest{
#' ##
#' X1<-matrix(rnorm(1000*2),1000,2)
#' simulated_data<-simfdata(n=1000,beta=c(0.5,0.5),fvar=0.5,
#' X=X1)
#' model1<-cphGM(formula=Surv(time,status)~X1+X2,
#' fterm<-c('gamma','id'),Time="time",status="status",
#' id="id",data=simulated_data,bhdist='weibull')
#' model1
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra K. Vishwakarma
#' @seealso \link{dscore},\link{simfdata}
#'
cphGM<-function(formula,fterm,Time,status,id,data,bhdist,method='L-BFGS-B',maxit=200){
  cl<-match.call()
  arg_checks<-as.list(match.call())[-1]
  arg_checks$formula<-deparse(formula)
  fterm_dist<-fterm[1]
  fterm_cluster<-fterm[2]
  c_id<-data[,fterm_cluster]
  # Check for null values
  if(!bhdist%in%c('weibull','exponential','gompertz')){
    stop('Available distribution are weibull,exponential,gompertz')
  }
  if(!fterm_dist%in%c('gamma')){
    stop('available frailty distribution is gamma')
  }
  if(is.null(data)) stop("data can't be NULL")
  if(is.null(status)) stop("status cannot be null")
  if(is.null(Time)) stop("Time cannot be null")
  if(is.null(id)) stop("id cannot be null")
  # Extract parts of the formula
  formula1<-deparse(formula)
  survmatch<-gregexpr("Surv\\([^()]+\\)",formula1)
  surv.part<-regmatches(formula1, survmatch)
  surv.statustime<-strsplit(gsub("Surv\\(([^()]+)\\)","\\1",surv.part),",")[[1]]
  surv.time<-gsub("\\s+","",surv.statustime[1])
  surv.status<-gsub("\\s+","", surv.statustime[2])
  surv.x<-strsplit(formula1,"\\~")[[1]][2]
  surv_x_var<-trimws(unlist(strsplit(surv.x,"\\+")))
  select_var<-apply(data[surv_x_var],2,function(x){is.numeric(x)})
  X_mean<-c()
  for(i in 1:length(select_var)){
    if(select_var[i]){
      data[,surv_x_var[i]]<-data[,surv_x_var[i]]-mean(data[,surv_x_var[i]])
      X_mean[i]<-mean(data[,surv_x_var[i]])
    }else{
      data[,surv_x_var[i]]<-data[,surv_x_var[i]]
      X_mean[i]<-NULL
    }
  }


  d_time<-data[,surv.time]
  unique_time<-sort(unique(d_time))
  d_status<-data[,surv.status]
  #X_cov<-model.matrix(as.formula(paste0("~",surv.x,"-1")),data)
  X_cov<-model.matrix(as.formula(paste0("~",surv.x)),data)
  X_cov<-X_cov[,-1]

  id<-data[,id]
  n_beta<-ncol(X_cov)


  #X_mean<-apply(X_cov,2,function(x){if(is.numeric(x)){mean(x)}else{NULL}})
  #X_cov<-apply(X_cov,2,function(x){if(is.numeric(x)){x<-scale(x,center=TRUE,scale=FALSE)}else{x<-x}})
  #

  # Define the partial likelihood function
  if(bhdist=='weibull'& fterm_dist=='gamma'){
    coxLikelihood<- function(param,Time,status,X,id){
      n<-nrow(X)
      bt<-param[1:ncol(X)]
      sig2<-exp(param[ncol(X)+1])
      p1<-exp(param[ncol(X)+2])
      p2<-exp(param[ncol(X)+3])
      Xbeta<-X%*%bt
      bh<-p1*p2*Time^(p1-1)
      bH<-p2*Time^p1
      log_likelihood<-0
      for(i in 1:n){
        log_likelihood_i<-status[i]*(log(bh[i])+Xbeta[i]-log(1+sig2*bH[i]*exp(Xbeta[i]))) -
          log(1+sig2*bH[i]*exp(Xbeta[i]))/sig2
        log_likelihood<-log_likelihood+log_likelihood_i
      }
      return(-log_likelihood)
    }
    #Initial parameter estimates
    param_init<-c(rep(0,n_beta),0.1,0.1,0.1)
    # Optimization to find the MLEs
    optResult<-optim(
      par=param_init,
      fn=coxLikelihood,
      Time=d_time,
      status=d_status,
      X=X_cov,
      id=id,
      method=method,
      control=list(maxit=maxit),
      hessian=TRUE
    )
    est_param<-optResult$par
    est_beta<-est_param[1:n_beta]
    est_theta<-exp(est_param[n_beta+1])
    est_rho<-exp(est_param[n_beta+2])
    est_lambda<-exp(est_param[n_beta+3])
    f_values<-sapply(unique(c_id), function(clusterID) {
      clstINDX<-which(c_id==clusterID)
      sum_status<-sum(d_status[clstINDX])
      sum_H<-sum(((est_lambda)*d_time[clstINDX]^est_rho)*exp(as.matrix(X_cov[clstINDX,])%*%est_beta))
      frailty_est<-(sum_status+(1/est_theta))/(sum_H+(1/est_theta))
      return(frailty_est)
    })
  }else if(bhdist=='exponential'& fterm_dist=='gamma'){

    coxLikelihood<- function(param,Time,status,X,id){
      n<-nrow(X)
      bt<-param[1:ncol(X)]
      sig2<-exp(param[ncol(X)+1])
      p1<-exp(param[ncol(X)+2])
      #p2<-exp(param[ncol(X)+3])
      Xbeta<-X%*%bt
      bh<-rep(p1,length(Time))
      bH<-p1*Time
      log_likelihood<-0
      for(i in 1:n){
        log_likelihood_i<-status[i]*(log(bh[i])+Xbeta[i]-log(1+sig2*bH[i]*exp(Xbeta[i]))) -
          log(1+sig2*bH[i]*exp(Xbeta[i]))/sig2
        log_likelihood<-log_likelihood+log_likelihood_i
      }
      return(-log_likelihood)
    }
    #Initial parameter estimates
    param_init<-c(rep(1,n_beta),0.1,0.1)
    # Optimization to find the MLEs
    optResult<-optim(
      par=param_init,
      fn=coxLikelihood,
      Time=d_time,
      status=d_status,
      X=X_cov,
      id=id,
      method=method,
      hessian=TRUE,
      control=list(maxit=maxit)
    )
    est_param<-optResult$par
    est_beta<-est_param[1:n_beta]
    est_theta<-exp(est_param[n_beta+1])#FRAILTY VARIANCE
    est_lambda<-exp(est_param[n_beta+2])  #LAMBDA: parameter for exponential distribution
    #est_lambda<-exp(est_param[n_beta+3])
    f_values<-sapply(unique(c_id), function(clusterID) {
      clstINDX<-which(c_id==clusterID)
      sum_status<-sum(d_status[clstINDX])
      sum_H<-sum((est_lambda*d_time[clstINDX])*exp(as.matrix(X_cov[clstINDX,])%*%est_beta))
      frailty_est<-(sum_status+(1/est_theta))/(sum_H+(1/est_theta))
      return(frailty_est)
    })

  }else if(bhdist=='gompertz'&fterm_dist=='gamma'){

    coxLikelihood<- function(param,Time,status,X,id){
      n<-nrow(X)
      bt<-param[1:ncol(X)]
      sig2<-exp(param[ncol(X)+1])#frailty variance
      p1<-exp(param[ncol(X)+2])#lambda
      p2<-exp(param[ncol(X)+3])#gamma
      Xbeta<-X%*%bt
      bh<-p1*exp(p2*Time)
      bH<-(p1/p2)*(exp(p2*Time)-1)
      log_likelihood<-0
      for(i in 1:n){
        log_likelihood_i<-status[i]*(log(bh[i])+Xbeta[i]-log(1+sig2*bH[i]*exp(Xbeta[i]))) -
          log(1+sig2*bH[i]*exp(Xbeta[i]))/sig2
        log_likelihood<-log_likelihood+log_likelihood_i
      }
      return(-log_likelihood)
    }
    #Initial parameter estimates
    param_init<-c(rep(1,n_beta),0.1,0.1,0.1)
    # Optimization to find the MLEs
    optResult<-optim(
      par=param_init,
      fn=coxLikelihood,
      Time=d_time,
      status=d_status,
      X=X_cov,
      id=id,
      method=method,
      control=list(maxit=maxit),
      hessian=TRUE
    )
    est_param<-optResult$par
    est_beta<-est_param[1:n_beta]
    est_theta<-exp(est_param[n_beta+1])
    est_lambda<-exp(est_param[n_beta+2])#lambda
    est_gamma<-exp(est_param[n_beta+3])#gamma
    f_values<-sapply(unique(c_id), function(clusterID) {
      clstINDX<-which(c_id==clusterID)
      sum_status<-sum(d_status[clstINDX])
      sum_H<-sum(((est_lambda/est_gamma)*exp(est_gamma*d_time[clstINDX]-1))*exp(as.matrix(X_cov[clstINDX,])%*%est_beta))
      frailty_est<-(sum_status+(1/est_theta))/(sum_H+(1/est_theta))
      return(frailty_est)
    })

  }

  result<-list()
  result$Call<-cl
  result$est_param<-est_param
  result$arguments<-arg_checks

  if(bhdist=='weibull'){
    result$est_beta<-est_beta
    result$est_theta<-est_theta
    result$est_rho<-est_rho
    result$est_lambda<-est_lambda
    b_haz<-est_lambda*unique_time^est_rho
    result$b_haz<-b_haz
  }else if(bhdist=='exponential'){
    result$est_beta<-est_beta
    result$est_theta<-est_theta
    result$est_lambda<-est_lambda
    b_haz<-est_lambda*unique_time
    result$b_haz<-b_haz
  }else if(bhdist=='gompertz'){
    result$est_beta<-est_beta
    result$est_theta<-est_theta
    result$est_lambda<-est_lambda
    result$est_gamma<-est_gamma
    b_haz<-(est_lambda/est_gamma)*exp((est_gamma*unique_time)-1)
    result$b_haz<-b_haz
  }

  result$frailty_values<-f_values
  result$value<-optResult$value
  result$convg<-optResult$convergence
  result$Hessian<-optResult$hessian
  result$SE<-diag(solve(optResult$hessian))
  result$b_hazard<-data.frame(Time=unique_time,C_Hazard=b_haz)
  #result$lp_pred<-X_cov%*%est_beta+log(rep(f_values,times=data.frame(table(c_id))$Freq))
  result$lp_pred<-X_cov%*%est_beta
  result$fterm_cluster<-fterm_cluster
  result$time<-Time
  result$status<-status
  result$id<-id
  result$cov<-X_cov
  result$X_mean<-X_mean
  class(result)<-'cphGM'
  return(result)
}

utils::globalVariables(c("model.matrix","as.formula","optim"))
