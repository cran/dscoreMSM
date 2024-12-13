#' @title Survival probability plot
#' @description it gives plot with fitted survival curve obtained from two different coxPH model fitted before and after SPSM
#' @param model1 coxPH fitted model object (before SPSM)
#' @param model2 coxPH fitted model object (after SPSM)
#' @param data1 multistate data used in model1
#' @param data2 multistate data used in model2
#' @param n_trans number of transition
#' @param id particular id from the dataset
#' @return plot for survival curve of a particular id obtained from both the model
#' @import ggplot2 mstate
#' @export
#' @examples
#' \donttest{
#' ##
#' library(mstate)
#' data(EBMTdata)
#' data(EBMTupdate)
#' tmat<-transMat(x=list(c(2,3),c(3),c()),names=c("Tx","Rec","Death"))
#' covs<-c("dissub","age","drmatch","tcd","prtime","x1","x2","x3","x4")
#' msbmt<-msprep(time=c(NA,"prtime","rfstime"),status=c(NA,"prstat","rfsstat"),
#'              data=EBMTdata,trans=tmat,keep=covs)
#' msbmt1<-msprep(time=c(NA,"prtime","rfstime"),status=c(NA,"prstat","rfsstat"),
#'               data=EBMTupdate,trans=tmat,keep=covs)
#' msph3<-coxph(Surv(time,status)~dissub+age+drmatch+tcd+
#'              frailty(id,distribution='gamma'),data=msbmt[msbmt$trans==3,])
#' msph33<-coxph(Surv(Tstart,Tstop,status)~dissub+age +drmatch+ tcd+
#'               frailty(id,distribution='gamma'),data=msbmt1[msbmt1$trans==3,])
#' ggplot_surv(model1=msph3,model2=msph33,data1=msbmt,
#'            data2=msbmt1,n_trans=3,id=1)
#' #####
#' # plot1<-ggplot_surv(model1=msph3,model2=msph33,data1=msbmt,data2=msbmt1,
#  # n_trans=3,id=1,folder_path="C:/Users/.....")
#' # ggsave("plot1.jpg",path="C:/Users/.....")
#' #####
#' ##
#' }
#' @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#' @seealso \link{dscore}, \link{simfdata}, \link{cphGM}
ggplot_surv<-function(model1,model2,data1,data2,n_trans,id){
  #msph1<-model1
  #msph11<-model2
  # Extract frailty terms for both models
  frailty_terms<-model1$frail
  frailty_terms_ps<-model2$frail
  # Get the frailty term for the specific individual (based on ID)
  frailty_ind<-frailty_terms[id]
  frailty_ind_ps<-frailty_terms_ps[id]
  # Calculate the baseline survival curve for both models
  base_surv<-survfit(model1)
  base_surv_ps<-survfit(model2)
  # Calculate the linear predictor for the individual, adding the frailty term
  lp_ind<-predict(model1,newdata=data1[data1$trans==n_trans&data1$id==id,],
                  type="lp",se.fit = FALSE)+frailty_ind
  lp_ind_ps<-predict(model2,newdata=data2[data2$trans==n_trans&data2$id==id,],
                     type="lp",se.fit = FALSE)+frailty_ind_ps
  # Compute the individual survival curves
  surv_ind<-base_surv$surv^exp(lp_ind)
  surv_ind_ps<-base_surv_ps$surv^exp(lp_ind_ps)
  # Extract time points for plotting
  time_points1<-base_surv$time
  time_points2<-base_surv_ps$time
  #pdf_filename<-paste0(folder_path,"/Survival_Plot_",n_trans,"_",id,".pdf")
  #pdf(file=pdf_filename)
  # Plot the individual survival curves
  df1 <- data.frame(
    Time = time_points1,
    Survival_Probability = surv_ind,
    Group = "No PSM"
  )

  df2 <- data.frame(
    Time = time_points2,
    Survival_Probability = surv_ind_ps,
    Group = "With PSM"
  )

  # Combine the data frames
  df <- rbind(df1, df2)

  # Create the ggplot
  p <- ggplot(df, aes(x = Time, y = Survival_Probability, color = Group, linetype = Group)) +
    geom_line(size = 1) +
    labs(
      title = paste("Individual Survival Curve (ID =", id, ")"),
      x = "Time",
      y = "Survival Probability"
    ) +
    scale_color_manual(values = c("blue", "red")) +
    scale_linetype_manual(values = c(1, 2)) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),legend.title = element_blank()) +
    ylim(0, 1)

  return(p)

}

utils::globalVariables(c("survfit","Time","Survival_Probability","Group"))
