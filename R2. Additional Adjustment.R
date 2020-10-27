library(dlnm)
library(mgcv)
library(mvmeta)

temp_bs <- O3_bs <- PM10_bs <- all.vcov.ts2.add <-vector("list",7)
all.coef.ts2.add <-matrix(NA,nrow=7,ncol=5)

cutoff_v <- c(seq(0.80,0.98,by=0.02)) ### Cut-off points for dummy variables
add.sens.all <- matrix(NA,nrow=length(cutoff_v),ncol=3)

for (iii in seq(length(cutoff_v))) {
  cutoff_perc <- cutoff_v[iii]
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
 
    
    ###Create additional dummy variables indicating time when these two model predict mortality differently
    ind.c<-c(rep(NA,365),as.numeric(abs((predict(all.fit.ts1))-(predict(all.fit.ts2)))
                        >quantile(abs((predict(all.fit.ts1))-(predict(all.fit.ts2))),cutoff_perc)
                                    
    ))
    ind.c[which(ind.c==1)]<-seq(length(ind.c[which(ind.c==1)]))
    
    
    ###Additional sensitivity analysis model
    all.fit.ts2.add<-glm(mt_b_total_all~PM10_bs[[aa]]+O3_bs[[aa]]+temp_bs[[aa]]+ns(rh_mean,3)+as.factor(weekday_2)+as.factor(HOLIDAY)
                                        +ns(week_to,3)+ns(month_to,3)+ns(DAYOFYEAR,20)+as.factor(year)+as.factor(inf_visit)+as.factor(ind.c),data=fit_dat,family=quasipoisson())
    
    
    all.coef.ts2.add[aa,1:5] <- coef(all.fit.ts2.add)[2:(5+1)]
    all.vcov.ts2.add[[aa]] <- vcov(all.fit.ts2.add)[2:(5+1),2:(5+1)]
    
    }
  
  ###POOLING
  all.mvmeta.ts2.add<-mvmeta(all.coef.ts2.add, all.vcov.ts2.add, method="reml", control=list(maxiter=1000))
  

  ###GET Results

  add.sens.all[iii,]<- (c(crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2.add),vcov=vcov(all.mvmeta.ts2.add),cen=0,model.link="log",lag=c(0,365),at=10)$allRRfit,
                          crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2.add),vcov=vcov(all.mvmeta.ts2.add),cen=0,model.link="log",lag=c(0,365),at=10)$allRRlow,
                          crosspred(PM10_bs[[1]],coef=coef(all.mvmeta.ts2.add),vcov=vcov(all.mvmeta.ts2.add),cen=0,model.link="log",lag=c(0,365),at=10)$allRRhigh)-1)*100
}  


###Display Results
par(mfrow=c(1,1),mar=c(4,3,0,0),mgp=c(2,1,0),oma=c(0.5,0.5,2,0.5))
plot(seq(10),rep(0,10),type="n",xlab="",ylab="Percent increase",xaxt="n",yaxt="n",ylim=c(-4,12),bty="n")
abline(h=0)
axis(2,seq(-4,12,by=4))
arrows(seq(10),add.sens.all[1:10,2],seq(10),add.sens.all[1:10,3],code=0)
points(seq(10),add.sens.all[1:10,1],pch=16,col="black")
axis(1,seq(10),seq(80,98,by=2),cex.axis=0.8,line=0)


