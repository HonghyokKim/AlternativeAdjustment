library(dlnm)
library(mgcv)
library(mvmeta)

temp_bs <- O3_bs <- PM10_bs <- all.vcov.ts1 <- all.vcov.ts2 <-  vector("list",7)

all.coef.ts1 <- all.coef.ts2 <- matrix(NA,nrow=7,ncol=5)

for ( aa in seq(7) ) {
    fit_dat<-data[[aa]]
    
    temp_bs[[aa]] <- crossbasis(data[[aa]]$ta_mean, lag=c(0,14), argvar=list(fun="bs",degree=2,
                                                                             knots=quantile(data[[aa]]$ta_mean,c(0.1,0.5,0.9))),
                                arglag=list(knots=logknots(21,fun="ns",4),fun="ns"))
    O3_bs[[aa]] <- crossbasis(data[[aa]]$o3_1000, lag=c(0,45), argvar=list(fun="lin"),arglag=list(knots=logknots(14,fun="ns",3),fun="ns"))
    PM10_bs[[aa]] <- crossbasis(data[[aa]]$PM10_c, lag=c(0,365), argvar=list(fun="lin"),arglag=list(knots=logknots(14,fun="ns",3),fun="ns"))
    
    ####NCS(t, 10df/year)
    all.fit.ts1<-glm(mt_b_total_all~PM10_bs[[aa]]+O3_bs[[aa]]+temp_bs[[aa]]+ns(rh_mean,3)+as.factor(weekday_2)+as.factor(HOLIDAY)
                                        +ns(TIME,10*8)+as.factor(inf_visit),data=fit_dat,family=quasipoisson())
    
    ####NCS(week,3df/year)+NCS(month,3df/year)+NCS(doy,20df)+I(year)
    all.fit.ts2<-glm(mt_b_total_all~PM10_bs[[aa]]+O3_bs[[aa]]+temp_bs[[aa]]+ns(rh_mean,3)+as.factor(weekday_2)+as.factor(HOLIDAY)
                                        +ns(week_to,3)+ns(month_to,3)+ns(DAYOFYEAR,20)+as.factor(year)+as.factor(inf_visit),data=fit_dat,family=quasipoisson())
 
    
    all.coef.ts1[aa,1:5] <- coef(all.fit.ts1)[2:(5+1)]
    all.vcov.ts1[[aa]] <- vcov(all.fit.ts1)[2:(5+1),2:(5+1)]
    all.coef.ts2[aa,1:5] <- coef(all.fit.ts2)[2:(5+1)]
    all.vcov.ts2[[aa]] <- vcov(all.fit.ts2)[2:(5+1),2:(5+1)]
    
  }
  
###POOLING
all.mvmeta.ts1<-mvmeta(all.coef.ts1, all.vcov.ts1 , method="reml", control=list(maxiter=1000))
all.mvmeta.ts2<-mvmeta(all.coef.ts2, all.vcov.ts2 , method="reml", control=list(maxiter=1000))
  


###GET Cumulative relative risk by exposure duration
all.CRR.mat.ts1 <- matrix(NA,nrow=366,ncol=3)
all.CRR.mat.ts2 <- matrix(NA,nrow=366,ncol=3)

for (aaa in seq(0,365)) {
  all.CRR.mat.ts1[aaa+1,]<-c(crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts1),vcov=vcov(all.mvmeta.ts1),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRfit,
                                     crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts1),vcov=vcov(all.mvmeta.ts1),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRlow,
                                     crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts1),vcov=vcov(all.mvmeta.ts1),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRhigh)
  all.CRR.mat.ts2[aaa+1,]<-c(crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2),vcov=vcov(all.mvmeta.ts2),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRfit,
                                     crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2),vcov=vcov(all.mvmeta.ts2),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRlow,
                                     crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2),vcov=vcov(all.mvmeta.ts2),cen=0,model.link="log",lag=c(0,aaa),at=10)$allRRhigh)
}


#### Percent increase by exposure duration ####
all.PI.mat.ts2 <- (all.CRR.mat.ts2-1)*100

par(mfrow=c(1,1),mgp=c(2,1,0),mar=c(4,4,1,1))
plot(x=seq(0,365),y=all.PI.mat.ts2[,1],type="l",xlab="Lag",ylab="Percent increase",bty="n",ylim=c(-10,20),xaxt="n")
abline(h=0,col="grey")
axis(1,at=seq(0,400,by=50))
polygon(c(rev(seq(0,365)), seq(0,365)), c(rev(all.PI.mat.ts2[,2]), all.PI.mat.ts2[,3]), col = rgb(0,0.8,1,0.5,maxColorValue = 1), border = NA)
lines(x=seq(0,365),y=all.PI.mat.ts2[,1],lty=1,col = rgb(0,0,0,0.7,maxColorValue = 1),lwd=3)

### RR by lag
par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(4,4,2,1))
plot(log="y",crossreduce(PM10_bs[[1]],coef=coef(all.mvmeta.ts1),vcov=vcov(all.mvmeta.ts1),type="var",cen=0,value=10,model.link="log"),
     ci.arg=list(col=rgb(0,0.8,1,0.5,maxColorValue = 1)),col =rgb(0,0,0,0.7,maxColorValue = 1),lwd=3,ylab="Relative Risk",ylim=c(0.998,1.004))
title("NCS(t,10/year)") ### overfitting

plot(log="y",crossreduce(PM10_bs[[1]],coef=coef(all.mvmeta.ts2),vcov=vcov(all.mvmeta.ts2),type="var",cen=0,value=10,model.link="log"),
     ci.arg=list(col=rgb(0,0.8,1,0.5,maxColorValue = 1)),col =rgb(0,0,0,0.7,maxColorValue = 1),lwd=3,ylab="Relative Risk",ylim=c(0.998,1.004))

title("Proposed") ### not overfitting



