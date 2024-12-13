#' print function for cphGM
#' @description
#' S3 print method for class 'cphGM'
#' @param x object
#' @param ... others
#' @return prints table containing various parameter estimates,
#'         SE, P-value.
#' @examples
#' \donttest{
#' ##
#' n1<-1000
#' p1<-2
#' X1<-matrix(rnorm(n1*p1),n1,p1)
#' simulated_data<-simfdata(n=1000,beta=c(0.5,0.5),fvar=0.5,X=X1)
#' model1<-cphGM(formula=Surv(time,status)~X1+X2,
#' fterm=c('gamma','id'),Time="time",status="status",
#' id="id",data=simulated_data,bhdist='weibull')
#' print(model1)
#' ##
#' }
#' @rdname print.cphGM
#' @method print cphGM
#' @export
print.cphGM<-function(x,...){
  x<-x
  digits<-3
  if(!inherits(x,'cphGM'))
    stop("\n Not a 'cphGM' object.\n")

  cat("\n coxPH model with parametric baseline hazard and gamma frailty")
  #cat("\n Call: \n",paste(deparse(x$Call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("\n Model formula:" ,x$arguments$formula ,'\n')
  cat("\n Baseline distribution:",x$arguments$bhdist,'\n')
  cat("\n Frailty distribution:", 'Gamma', "\n")
  cat("\n Frailty variable:",x$fterm_cluster,"\n")

  cat('\n Estimates: \n')
  vname<-colnames(x$cov)


  #lstan<-x$pseudoMod$stan_summary[starts_with('Long1',vars =rownames(x$pseudoMod$stan_summary)),c(1,3,4,10)]
  #lstan<-cbind(lstan,Zvalue=lstan[,1]/lstan[,2])
  f1<-function(x){
    return(if((x[1]>0))2*(1-pnorm((x[1]/x[2]),mean=0,sd=1))else 2*(pnorm((x[1]/x[2]),mean=0,sd=1)))
  }

  if(x$Call$bhdist=='weibull'){
  est<-c(x$est_beta,x$est_theta,x$est_rho,x$est_lambda)
  est_se<-sqrt(x$SE)
  zvalue<-est/est_se
  est_table<-data.frame(est,est_se)
  Pvalue<-apply(est_table,1,f1)
  est_table<-cbind(est_table,zvalue,Pvalue)
  attr(est_table,'names')<-c('Est','StDev','Z-value','P-value')
  attr(est_table,'row.names')<-c(vname,'Theta','Rho','Lambda')
  }else if(x$Call$bhdist=='exponential'){
    est<-c(x$est_beta,x$est_theta,x$est_lambda)
    est_se<-sqrt(x$SE)
    zvalue<-est/est_se
    est_table<-data.frame(est,est_se)
    Pvalue<-apply(est_table,1,f1)
    est_table<-cbind(est_table,zvalue,Pvalue)
    attr(est_table,'names')<-c('Est','StDev','Z-value','P-value')
    attr(est_table,'row.names')<-c(vname,'Theta','Lambda')
  }else{
    est<-c(x$est_beta,x$est_theta,x$est_lambda,x$est_gamma)
    est_se<-sqrt(x$SE)
    zvalue<-est/est_se
    est_table<-data.frame(est,est_se)
    Pvalue<-apply(est_table,1,f1)
    est_table<-cbind(est_table,zvalue,Pvalue)
    attr(est_table,'names')<-c('Est','StDev','Z-value','P-value')
    attr(est_table,'row.names')<-c(vname,'Theta','Lambda','Gamma')
  }

  est_table<-round(est_table,digits=digits)

  print(est_table)

  invisible(x)
}

