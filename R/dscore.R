#' @title Survival Proximity Score matching for MSM
#' @description function for survival proximity score matching in multistate model with three state.
#' @param status status column name in the survival data
#' @param data survival data
#' @param prob threshold probability
#' @param m starting column number
#' @param n ending column number
#' @param method distance metric name e.g. "euclidean","minkowski","canberra"
#' @return list with newdataset updated using dscore
#' @import stats
#' @export
#' @references
#' Vishwakarma, G. K., Bhattacherjee, A., Rajbongshi, B. K., & Tripathy, A. (2024). Censored imputation of time to event outcome through survival proximity score method. \emph{Journal of Computational and Applied Mathematics}, 116103;
#'
#' Bhattacharjee, A., Vishwakarma, G. K., Tripathy, A., & Rajbongshi, B. K. (2024). Competing risk multistate censored data modeling by propensity score matching method. \emph{Scientific Reports}, 14(1), 4368.
#' @examples
#' \donttest{
#' ##s
#' data(simulated_data)
#' udata<-dscore(status="status",data=simulated_data,prob=0.65,m=4,n=7)
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra K. Vishwakarma
#' @seealso \link{cphGM},\link{simfdata}
dscore<-function(status,data,prob,m,n,method="euclidean"){
  # dataset must contain two columns named status,time
  mdata<-data
  mdata[is.na(mdata)]<-0
  k<-which(colnames(mdata)==status)
  mdata1<-mdata[,m:n]
  mdata2<-cbind(status=mdata$status,time=mdata$time,mdata[,m:n])
  d11<-subset(mdata2,mdata2$status==1)
  d00<-subset(mdata2,mdata2$status==0)
  coefficient<-c()
  for(i in 1:ncol(mdata2)-2){
    fit<-survmc(m=i+2,n=i+2,Time="time",Event="status",chains=2,adapt=100,iter=50,data=mdata2)
    if(fit$`Posterior Means`&fit$`2.5%`&fit$`97.5%`>1){
      coefficient[i]<-fit$`Posterior Means`
    }else if(fit$`Posterior Means`&fit$`2.5%`&fit$`97.5%`<1){
      coefficient[i]<-fit$`Posterior Means`
    } else {
      coefficient[i]<-0
    }
  }
  #parm<-summary(model1)$coefficients[,1]
  d1<-d11[,-which(colnames(d11)=="status")]
  t<-which(colnames(d1)=="time")
  d1<-d1[,-t]
  d0<-d00[,-which(colnames(d11)=="status")]
  xmin<-c()
  k<-length(coefficient)

  xmin<-function(x,data){
    xmin<-c()
    for(i in 1:length(x)){
      if(x[i]>0){
        xmin[i]<-min(data[,i])
      }else{xmin[i]<-max(data[,i])}
    }
    xmin
  }
  xmin1<-xmin(x=coefficient,data=d1)


  Dscore<-c()
  for( i in 1:nrow(d0)){
    if(anyNA(xmin1)==T){
      k<-which(is.na(xmin1)==T)
    }
    d<-rbind(d0[i,][-k],xmin1[-k])
    Dscore[i]<-dist(d,method=method)
  }
  rslt<-list()
  rslt$dscore<-Dscore
  rslt$olddeath<-nrow(d11)
  risk<-ecdf(Dscore)
  proximity<-1-risk(Dscore)
  nd00<-cbind(d00,proximity)
  nstatus<-c()
  for(i in 1:nrow(nd00)){
    if(nd00$proximity[i]>prob){
      nstatus[i]<-1
    }else{
      nstatus[i]<-0
    }
  }
  d11<-subset(mdata,mdata$status==1)
  d00<-subset(mdata,mdata$status==0)
  d00$status<-as.factor(nstatus)
  newdata<-rbind(d11,d00)
  newdata<-newdata[order(newdata$id),]
  rslt$updateddeath<-nrow(subset(newdata,newdata$status=='1'))
  rslt$newdata<-newdata
  #-------------------------------------------------------------------------
  rslt
}

utils::globalVariables(c("dist","ecdf"))
