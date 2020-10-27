

simtotal<-10000
numlags<-730
ptype<-"Type1"
concurvity<-1

X_IRR <- vector("list",5)
X_CRR <- vector("list",5)
X_CRR_95CI <- vector("list",5)
X_CRR_COVERAGE <- vector("list",5)
for (aa in seq(5)) {
    X_IRR[[aa]]<- matrix(NA,nrow=simtotal,ncol=numlags+1)
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

    X_bs <- crossbasis(X, lag=c(0,730),argvar=list(fun="lin"),arglag=list(fun="ns",knots=logknots(730,fun="ns",10) ) )
    
    #####FINAL DATASET OBS = 3653.  BUT the front 21 obs of newcount is missing. ######
    for (ss in seq(4)) {
      
      if (ss == 1) { ## SEASONALITY & LONG-TERM TIME-TREND ARE NOT ADJUSTED
        model_X_fit<-glm(newcount~X_bs[(365*2+1):(4383),]+tempbs[(365*2+1):(4383),]+
                           ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit), data=newdata, family=poisson())
      } 
      
    
      if (ss == 2) { ## PROPOSED METHOD
        model_X_fit<-glm(newcount~X_bs[(365*2+1):(4383),]+tempbs[(365*2+1):(4383),]+
                           ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit)+
                           week_o_1+week_o_2+week_o_3+week_o_4+week_o_5+
                           month_o_1+month_o_2+month_o_3+month_o_4+month_o_5+as.factor(year)+ns(DAYOFYEAR,10), data=newdata, family=poisson())
      }   
      
      if (ss == 3) { ## STANDARD METHOD NCS(t,4/year)
        model_X_fit<-glm(newcount~X_bs[(365*2+1):(4383),]+tempbs[(365*2+1):(4383),]+
                           ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit)+
                           ns(TIME,4*10)+as.factor(inf_visit), data=newdata, family=poisson())
      }   
      
      if (ss == 4) { ## STANDARD METHOD NCS(t,10/year)
        model_X_fit<-glm(newcount~X_bs[(365*2+1):(4383),]+tempbs[(365*2+1):(4383),]+
                           ns(rh_mean,4)+as.factor(HOLIDAY)+as.factor(weekday_2)+as.factor(inf_visit)+
                           ns(TIME,10*10)+as.factor(inf_visit), data=newdata, family=poisson())
      }   
      
    
      
      X_IRR[[ss]][nsim,1:(numlags+1)] <- crossreduce(X_bs,model_X_fit, type="var", value=1, cen=0)$RRfit
      X_CRR[[ss]][nsim,] <- crossreduce(X_bs,model_X_fit, at=1)$RRfit
      X_CRR_95CI[[ss]][nsim,1:2] <- c(crossreduce(X_bs,model_X_fit, at=1)$RRlow,crossreduce(X_bs,model_X_fit, at=1)$RRhigh)
      coverage <- findInterval(
        exp(sum(Beta[[ptype]][1:(numlags+1),1])),
        X_CRR_95CI[[ss]][nsim,1:2],
        left.open=TRUE,rightmost.close=TRUE)
      X_CRR_COVERAGE[[ss]][nsim,] <- ifelse(coverage==1,1,0)
    }  

    end.time <- Sys.time()
    print(paste(nsim,"run", round(nsim/simtotal*100,1),"% completed," ,round(end.time-start.time,3),"sec"))
}




####Summary####
CRR_summary <- matrix(NA,nrow=4,ncol=3)
colnames(CRR_summary) <- c("Bias","SD1000","Coverage")
rownames(CRR_summary) <- c("Noadj","PROPOSED","NCS(t,4/year)","NCS(t,10/year)")


for (ss in seq(4)) {
  CRR_summary[ss,1] <-(mean(log(X_CRR[[ss]]))-sum(Beta[[ptype]][1:(numlags+1),1]))/sum(Beta[[ptype]][1:(numlags+1),1])*100
  CRR_summary[ss,2] <-sd(log(X_CRR[[ss]]))*1000
  CRR_summary[ss,3] <-mean(X_CRR_COVERAGE[[ss]])
}

CRR_summary


CRR_TRACE <- matrix(NA,nrow=simtotal,ncol=1)
for (ii in seq(simtotal)) {
  
  CRR_TRACE[ii,1] <- mean(log(X_CRR[[2]])[1:ii])
}
par(mfrow=c(2,1))
plot((CRR_TRACE),xlab="Simulation #",ylab="Average of log-CRR",type="l")
abline(h=(sum(Beta[[ptype]][,1])),col="grey")

plot(colMeans(X_IRR[[2]], na.rm=T), type="n", ylim=c(0.9999,1.0003), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,730), xaxt="n")
lines(exp(Beta[[ptype]][,1]),col="grey", lwd=3, lty=1)

lines(exp(colMeans(log(X_IRR[[2]]),na.rm=T)), col="blue",lwd=1,lty=1)
abline(h=1)



