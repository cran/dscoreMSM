# @title Reciever Operating Curve
# @description this function provides roc plot for coxph model fitted before and after survival proximity score matching.
# @param trns transition number for the multistate model
# @param model1 fitted object from coxPH before SPSM
# @param model2 fitted object from coxPH after SPSM
# @param data1 dataset used for model1
# @param data2 dataset used for model2
# @param folder_path folder path where the user wants to save the plots
# @return returns roc plot for model1 and model2
# @import timeROC survival mstate
# @references
# Vishwakarma, G. K., Bhattacherjee, A., Rajbongshi, B. K., & Tripathy, A. (2024). Censored imputation of time to event outcome through survival proximity score method. \emph{Journal of Computational and Applied Mathematics}, 116103;
#
# Bhattacharjee, A., Vishwakarma, G. K., Tripathy, A., & Rajbongshi, B. K. (2024). Competing risk multistate censored data modeling by propensity score matching method. \emph{Scientific Reports}, 14(1), 4368.
# @examples
# \donttest{
# ##
# library(mstate)
# data(EBMTdata)
# data(EBMTupdate)
# tmat<-transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "Rec", "Death"))
# covs<-c("dissub", "age", "drmatch", "tcd", "prtime","x1","x2","x3","x4")
# msbmt<-msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
#              data = EBMTdata, trans = tmat, keep = covs)
# msbmt1<-msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"),
#               data = EBMTupdate, trans = tmat, keep = covs)
#
# msph3<-coxph(Surv(time,status)~dissub+age +drmatch+ tcd+
#              frailty(id,distribution = 'gamma'),data=msbmt[msbmt$trans==3,])
# msph33<-coxph(Surv(Tstart,Tstop,status)~dissub+age +drmatch+ tcd+
#               frailty(id,distribution = 'gamma'),data=msbmt1[msbmt1$trans==3,])
#
# #roc_plot(trns=3,model1=msph3,model2=msph33,
# #folder_path="C:/Users/bhrigu/OneDrive/Desktop/Paper5/docpresentation/euclidean",
# #data1=msbmt,data2=msbmt1)
# ##
#
# }
# @author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
# @seealso \link{dscore}, \link{simfdata}, \link{cphGM}
roc_plot <- function(trns, model1, model2, data1, data2, folder_path) {
  pred_model1 <- predict(model1, type = "lp")
  pred_model2 <- predict(model2, type = "lp")
  roc1 <- timeROC(T = data1[data1$trans == trns,]$time,
                  delta = data1[data1$trans == trns,]$status,
                  marker = pred_model1,
                  cause = 1,
                  times = quantile(data1[data1$trans == trns,]$time, probs = seq(0.2, 0.8, 0.1)))
  roc2 <- timeROC(T = data2[data2$trans == trns,]$time,
                  delta = data2[data2$trans == trns,]$status,
                  marker = pred_model2,
                  cause = 1,
                  times = quantile(data2[data2$trans == trns,]$time, probs = seq(0.2, 0.8, 0.1)))

  auc1 <- roc1$AUC
  auc2 <- roc2$AUC
  time_seq <- quantile(data2[data2$trans == trns,]$time, probs = seq(0.2, 0.8, 0.1))

  pdf_filename <- paste0(folder_path, "/ROC_Plots","_",trns,".pdf")
  pdf(file = pdf_filename)
  for (i in 1:length(time_seq)) {
    plot(roc1, time = time_seq[i], col = "blue", title = FALSE, legend = FALSE)
    plot(roc2, time = time_seq[i], add = TRUE, col = "red")
    legend("bottomright", legend = c("No PSM", "with PSM"), col = c("blue", "red"), lwd = 2)
    title(main = paste("Time=", round(time_seq[i],2), "AUC(", round(auc1[i], 2), ",", round(auc2[i], 2), ")"))
  }
  dev.off()
}

utils::globalVariables(c("pdf","legend","title","dev.off"))
