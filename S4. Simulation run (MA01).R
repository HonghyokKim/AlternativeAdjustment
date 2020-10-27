

simtotal<-5000
ptype<-"Type2"
concurvity<-1




X_CRR <- vector("list",17)
X_CRR_95CI <- vector("list",17)
X_CRR_COVERAGE <- vector("list",17)
for (aa in seq(17)) {
  X_CRR[[aa]] <- matrix(NA,nrow=simtotal,ncol=1)
  X_CRR_95CI[[aa]] <- matrix(NA,nrow=simtotal,ncol=2)
  X_CRR_COVERAGE[[aa]] <- matrix(NA,nrow=simtotal,ncol=1)
}


##############Simulation runs########
for (nsim in seq(simtotal)) {
  set.seed(nsim)
  start.time <- Sys.time()
  
  PredX <- PM10.pred
  PredY <- Y.pred
  
  Yset<-c(rep(NA,21),PredY)
  Xset <- PredX
  Xset_resid_sigma <- sigma(logPM10.fit)/concurvity
  
  X <- exp( rnorm(length(Xset), log(Xset), sd=Xset_resid_sigma) )
  X <- c(rep(NA,21),X)
  Xset_use <- c(rep(NA,21),Xset) 
  
  newdata <- sim_dat[(731):(4383),]
  newdata$X <- X[(731):(4383)]
  newdata$Xset <- Xset[(731):(4383)]
  
  X <- X-mean(X,na.rm=T)
  Xlag <- matrix(NA,nrow=length(X),ncol=T2length)
  for (j in seq(T2length)) {
    Xlag[,j] <- c(Lag(X,j-1))
  }
  XBeta <- Xlag%*%Beta[[ptype]][,1]
  
  Yset<-Yset[(731):(4383)]
  Xbetaset<-exp(XBeta[(731):(4383)])
  Lamda<-Yset*Xbetaset
  
  ####  obs = 4383 - 365 *2 (front) = 3653 ###
  newdata$newcount <-c(rep(NA,21),rpois(length(Lamda[c(-21:-1)]),Lamda[c(-21:-1)]))
  
  ##### obs = 4383 ###
  tempbs <- crossbasis(sim_dat$ta_mean, lag=c(0,21),argvar=list(fun="bs",degree=2,knots=tempknots),arglag=list(knots=logknots(21, fun="bs", df=7)))
  
  ##### obs = 4383, but the front 21 obs is missing ###
  
  newdata$X_ma01 <- runMean(newdata$X,0:1)
  
  X2_bs <- crossbasis(X, lag=c(2,730),argvar=list(fun="lin"),arglag=list(fun="ns",knots=seq(2,730,by=100) ) )
  
  
  #####FINAL DATASET OBS = 3653.  BUT the front 21 obs of newcount is missing. ######
  for (ss in seq(17)) {
    
    if (ss == 1) { ## SEASONALITY & LONG-TERM TIME-TREND ARE NOT ADJUSTED
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit), data=newdata, family=poisson())
    } 
    
    if (ss %in% c(2:13)) { ## STANDARD METHOD NCS(t,4-15/year)
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit)+
                         ns(TIME,(ss+2)*10)+as.factor(inf_visit), data=newdata, family=poisson())
    }   

    if (ss == 14) { ## STANDARD METHOD NCS(20/year)
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit)+
                         ns(TIME,20*10)+as.factor(inf_visit), data=newdata, family=poisson())
    }   
    
    if (ss == 15) { ## SEASONALITY & LONG-TERM TIME-TREND ARE NOT ADJUSTED BUT DL2-730 ADJUSTED
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+X2_bs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit), data=newdata, family=poisson())
    } 
    
    if (ss == 16) { ## NCS(10/year) + DL2-730 ADJUSTED
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+X2_bs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+ns(TIME,10*10)+as.factor(inf_visit), data=newdata, family=poisson())
    }
    if (ss == 17) { ## NCS(20/year) + DL2-730 ADJUSTED
      model_X_fit<-glm(newcount~X_ma01+tempbs[(365*2+1):(4383),]+X2_bs[(365*2+1):(4383),]+
                         ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+ns(TIME,20*10)+as.factor(inf_visit), data=newdata, family=poisson())
    }
    
    
    X_CRR[[ss]][nsim,] <- exp(coef(model_X_fit)[2])
    X_CRR_95CI[[ss]][nsim,1:2] <- exp(c(coef(model_X_fit)[2]-1.96*sqrt(vcov(model_X_fit)[2,2]),coef(model_X_fit)[2]+1.96*sqrt(vcov(model_X_fit)[2,2])))
    coverage <- findInterval(
      exp(sum(Beta[[ptype]][1:2,1])),
      X_CRR_95CI[[ss]][nsim,1:2],
      left.open=TRUE,rightmost.close=TRUE)
    X_CRR_COVERAGE[[ss]][nsim,] <- ifelse(coverage==1,1,0)
    
  }  
  
  end.time <- Sys.time()
  print(paste(nsim,"run", round(nsim/simtotal*100,1),"% completed," ,round(end.time-start.time,3),"sec"))
}




####Summary####
CRR_summary <- matrix(NA,nrow=17,ncol=3)
colnames(CRR_summary) <- c("Bias","SD","Coverage")
rownames(CRR_summary) <- c("Noadj",paste0("NCS(",seq(4,15),"/year)"),"NCS(20/year)","Noadj + DL2_730","NCS(10/year) + DL2_730","NCS(20/year) + DL2_730")


for (ss in seq(17)) {
  CRR_summary[ss,1] <-(mean(log(X_CRR[[ss]]))-sum(Beta[[ptype]][1:2,1]))/sum(Beta[[ptype]][1:2,1])*100
  CRR_summary[ss,2] <-sd(log(X_CRR[[ss]]))
  CRR_summary[ss,3] <-mean(X_CRR_COVERAGE[[ss]])
}

CRR_summary



