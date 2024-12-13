#@title survmc
#@description Bayesian Survival Model for high dimensional survival survival data
#@param m Starting column number from where variables of high dimensional data will get selected.
#@param n Ending column number till where variables of high dimensional data will get selected.
#@param Time Variable/Column name containing the information on duration of survival
#@param Event Variable/Column name containing the information of survival event
#@param chains Number of chains to perform
#@param adapt Number of adaptations to perform
#@param iter Number of iterations to perform
#@param data High dimensional data having survival duration and event.
#@return returns information on Bayesian survival model output
#@import rjags
#@references Rizopoulos, D., G. Papageorgiou, and P. Miranda Afonso. "JMbayes2: extended joint models for longitudinal and time-to-event data." R package version 0.2-4 (2022).
#@examples
#\donttest{
###
#1+2
#1+2
###
#}
#@author Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra Kumar Vishwakarma
#@seealso \link{dscore}
#' @import rjags
survmc <- function(m,n,Time,Event,chains,adapt,iter,data)
{
  if(Time!="OS"){
    names(data)[names(data) == Time] <- "OS"
  }
  if(Event!="Death"){
    names(data)[names(data) == Event] <- "Death"
  }
  data<-data[order(data$OS),]
  var1 <- colnames(data)
  nr <- nrow(data)
  data1 <- subset(data, Death == 1) #subsetting data with death status = 1
  u <- unique(data1$OS) #creating a vector with unique values of OS
  #adding a condition for censoring time vector to include the last censored patient when censoring = 0
  if ((data$Death[nrow(data)])==0){
    u1<-c(u,data$OS[nrow(data)])
  } else {
    u1 <- u
  }
  u2 <- sort(u1)
  u2[length(u2)]<-u2[length(u2)]+1E-3
  t.len<-(length(u2)-1)
  model_jags <- "
  data{
    # Set up data
  for(i in 1:N) {
    for(j in 1:T) {
    Y[i,j] <- step(obs.t[i] - t[j] + eps)
    dN[i, j] <- Y[i, j] * step(t[j + 1] - obs.t[i] - eps) * fail[i]
    }
  }
  }

  # Model
  model{
  for(i in 1:N){
    betax[i,1] <- 0
    for(k in 2:(p+1)){
      betax[i,k] <- betax[i,k-1] + beta[k-1]*x[i,k-1]
    }
  }
  for(j in 1:T) {
    for(i in 1:N) {
    dN[i, j] ~ dpois(Idt[i, j]) # Likelihood
    Idt[i, j] <- Y[i, j] * exp(betax[i,p+1]) * dL0[j] # Intensity
    }
    dL0[j] ~ dgamma(mu[j], c)
    mu[j] <- dL0.star[j] * c # prior mean hazard
  }
  c <- 0.001
  r <- 0.1
  for (j in 1 : T) {
    dL0.star[j] <- r * (t[j + 1] - t[j])
  }
  for(k in 1:p){
    beta[k] ~ dnorm(0.0,0.000001)
  }
  }"

  params <- c("beta","dL0")
  inits <-  function(){list( beta = rep(0,p), dL0 = rep(0.0001,bigt))}
  x2 <- rep(0,nrow(data))
  q <- matrix(nrow=0,ncol=5)
  s <- matrix(nrow=0,ncol=2)
  di <- matrix(nrow=0,ncol=1)
  for(i in m:n){
    x1 <- data[(1:nrow(data)),i]
    x = t(rbind(x1,x2))
    datafi <- list(x=x,obs.t=data$OS,t=u2,T=t.len,N=nrow(data),fail=data$Death,eps=1E-10,p=2)
    jags <- jags.model(textConnection(model_jags),
                       data = datafi,
                       n.chains = chains,
                       n.adapt = adapt)
    samps <- coda.samples(jags, params, n.iter=iter)
    s1 <- summary(samps)
    stats <- s1$statistics[1,c(1:2)]
    s <- rbind(s,stats)
    quan <- s1$quantiles[1,]
    q <- rbind(q,quan)
    d = dic.samples(jags, n.iter=iter)
    meandeviance <- round(sum(d$deviance),2)
    di <- rbind(di,meandeviance)
  }
  results <- cbind(s,q)
  expresults <- exp(results)
  Variables <- names(data)[m:n]
  expresults <- data.frame(Variables,expresults,di)
  colnames(expresults)<-c("Variable","Posterior Means","SD","2.5%","25%","50%","75%","97.5%","DIC")
  rownames(expresults) <- NULL
  return(expresults)
}

utils::globalVariables(c("Death","p","bigt","jags.model","coda.samples","dic.samples"))

