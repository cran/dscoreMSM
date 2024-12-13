## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## -----------------------------------------------------------------------------
library(dscoreMSM)
library(mstate)
library(rjags)
library(ggplot2)
library(survival)

## -----------------------------------------------------------------------------
n1<-1000
p1<-2
X1<-matrix(rnorm(n1*p1),n1,p1)
simulated_data<-simfdata(n=1000,beta=c(0.5,0.5),fvar=0.5,X=X1)

## -----------------------------------------------------------------------------
model1<-cphGM(formula=Surv(time,status)~X1+X2,
              fterm=c('gamma','id'),Time="time",status="status",
              id="id",data=simulated_data,bhdist='weibull')

## -----------------------------------------------------------------------------
data.frame("True value"=c(0.5,0.5,0.5,2,1),"Est value"=
             c(model1$est_beta,
               model1$est_theta,
               model1$est_rho,
               model1$est_lambda))

## -----------------------------------------------------------------------------
hist(model1$frailty_values,main="Histogram for z_i",xlab="estimated frailty")

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
msm_data<-function(){
tmat <- transMat(x = list(c(2, 3), c(3), 
                          c()), names = c("Tx", "Rec", "Death"))
covs<-data.frame(x1=rnorm(100,0,1),x2=rnorm(100,10,20),x3=rnorm(100,5,4),x4=rnorm(100,3,2))
covs<-as.matrix(covs)
ssim13<-simfdata(n=100,beta=c(0.1,0.1,0.5,0.5),fvar=0.5,X=covs)
ssim23<-simfdata(n=100,beta=c(0.2,0.02,0.5,0.04),fvar=0.5,X=covs)

stime13<-ssim13$time
stime23<-ssim23$time
#sstatus12<-ssim12$status
sstatus13<-ssim13$status
sstatus23<-ssim23$status
data1<-data.frame(id=1:100, stime13,stime23,sstatus13,sstatus23,covs)

for(i in 1:nrow(data1)){
  if(data1[i,]$sstatus13==0&data1[i,]$sstatus23==0){
    data1[i,]$stime23=data1[i,]$stime13
  }else if(data1[i,]$sstatus13==1&data1[i,]$sstatus23==0) {
    data1[i,]$stime23=data1[i,]$stime13+data1[i,]$stime23
  } else if(data1[i,]$sstatus13==1&data1[i,]$sstatus23==1){
    data1[i,]$stime23=data1[i,]$stime13+data1[i,]$stime23
  } else if(data1[i,]$sstatus13==0&data1[i,]$sstatus23==1){
    data1[i,]$stime13=data1[i,]$stime23
  }
}
udata1<-list()
udata1$id<-data1$id;udata1$status<-data1$sstatus23;udata1$time<-data1$stime23;udata1$x1<-data1$x1;udata1$x2<-data1$x2
udata1$x3<-data1$x3;udata1$x4<-data1$x4
udata1<-as.data.frame(udata1)
udata<-dscore(status="status",data=udata1,prob=0.65,m=4,n=7)
ndata<-udata$newdata
data2<-data1
data2$sstatus23<-as.numeric(as.character(ndata$status))
covs1 <- c("x1", "x2", "x3","x4")
msbmt <- msprep(time = c(NA, "stime13", "stime23"), status = c(NA,
                                                               "sstatus13", "sstatus23"), data =data1, trans = tmat, keep = covs1)

msbmt1 <- msprep(time = c(NA, "stime13", "stime23"), status = c(NA,
                                                                "sstatus13", "sstatus23"), data =data2, trans = tmat, keep = covs1)

msbmt <- expand.covs(msbmt, covs1, append = TRUE, longnames = FALSE)

msbmt1 <- expand.covs(msbmt1, covs1, append = TRUE, longnames = FALSE)
result<-list()
result$msbmt<-msbmt
result$msbmt1<-msbmt1
return(result)
}
sim_data<-msm_data()
example_data1<-sim_data[[1]]
example_data1_updated<-sim_data[[2]]
table(example_data1$status)
table(example_data1_updated$status)

## ----eval=FALSE---------------------------------------------------------------
# data("EBMTdata")
# data("EBMTupdate")
# table(EBMTdata$status)
# table(EBMTupdate$status)
# tmat<-transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "Rec", "Death"))
# covs<-c("dissub", "age", "drmatch", "tcd", "prtime","x1","x2","x3","x4")
# msbmt<-msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
#               data = EBMTdata, trans = tmat, keep = covs)
# msbmt1<-msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
#                data = EBMTupdate, trans = tmat, keep = covs)
# 
# msph3<-coxph(Surv(time,status)~dissub+age +drmatch+ tcd+
#                frailty(id,distribution = 'gamma'),data=msbmt[msbmt$trans==3,])
# msph33<-coxph(Surv(Tstart,Tstop,status)~dissub+age +drmatch+ tcd+
#                 frailty(id,distribution = 'gamma'),data=msbmt1[msbmt1$trans==3,])
# 
# ggplot_roc(trns=3,model1=msph3,model2=msph33,
#            folder_path=NULL,
#            data1=msbmt,data2=msbmt1,times=NULL)

## ----eval=FALSE---------------------------------------------------------------
# ggplot_surv(model1=msph3,model2=msph33,data1=msbmt,data2=msbmt1,n_trans=3,id=1)

