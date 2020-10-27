
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE,cex=0.8)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#####RR lags Figure

load("save simulation results and load here.Rdata")
LowConcurvity.Beta1_CRRsummary <- CRR_summary
LowConcurvity.Beta1_IRR <- X_IRR

load("save simulation results and load here.Rdata")
LowConcurvity.Beta2_CRRsummary <- CRR_summary
LowConcurvity.Beta2_IRR <- X_IRR

load("save simulation results and load here.Rdata")
HighConcurvity.Beta1_CRRsummary <- CRR_summary
HighConcurvity.Beta1_IRR <- X_IRR

load("save simulation results and load here.Rdata")
HighConcurvity.Beta2_CRRsummary <- CRR_summary
HighConcurvity.Beta2_IRR <- X_IRR





#jpeg("FILE FOLDER AND FILE NAME.jpg",
#     width=10, height=12, pointsize=15,res=900,units="in")
par(mfrow=c(4,2),mgp=c(2,1,0),mar=c(3.5,4,3,1),oma=c(5,0,0,0))
plot(colMeans(LowConcurvity.Beta1_IRR[[1]], na.rm=T), type="n", ylim=c(0.999,1.002), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[1]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(LowConcurvity.Beta1_IRR[[1]])*10,na.rm=T)), col="royalblue1",lwd=2,lty=2)
lines(exp(colMeans(log(LowConcurvity.Beta1_IRR[[2]])*10,na.rm=T)), col="darkblue",lwd=2,lty=3)
text(-140,1.00275,"A)",cex=1.2,xpd=TRUE)

plot(colMeans(LowConcurvity.Beta1_IRR[[1]], na.rm=T), type="n", ylim=c(0.990,1.010), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[1]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(LowConcurvity.Beta1_IRR[[3]])*10,na.rm=T)), col="lightskyblue1",lwd=2,lty=2)
lines(exp(colMeans(log(LowConcurvity.Beta1_IRR[[4]])*10,na.rm=T)), col="lightskyblue4",lwd=2,lty=3)
text(-140,1.015,"B)",cex=1.2,xpd=TRUE)


plot(colMeans(LowConcurvity.Beta2_IRR[[1]], na.rm=T), type="n", ylim=c(0.999,1.002), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[2]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(LowConcurvity.Beta2_IRR[[1]])*10,na.rm=T)), col="royalblue1",lwd=2,lty=2)
lines(exp(colMeans(log(LowConcurvity.Beta2_IRR[[2]])*10,na.rm=T)), col="darkblue",lwd=2,lty=3)
text(-140,1.00275,"C)",cex=1.2,xpd=TRUE)

plot(colMeans(LowConcurvity.Beta2_IRR[[1]], na.rm=T), type="n", ylim=c(0.990,1.010), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[2]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(LowConcurvity.Beta2_IRR[[3]])*10,na.rm=T)), col="lightskyblue1",lwd=2,lty=2)
lines(exp(colMeans(log(LowConcurvity.Beta2_IRR[[4]])*10,na.rm=T)), col="lightskyblue4",lwd=2,lty=3)
text(-140,1.015,"D)",cex=1.2,xpd=TRUE)

plot(colMeans(HighConcurvity.Beta1_IRR[[1]], na.rm=T), type="n", ylim=c(0.999,1.002), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[1]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(HighConcurvity.Beta1_IRR[[1]])*10,na.rm=T)), col="royalblue1",lwd=2,lty=2)
lines(exp(colMeans(log(HighConcurvity.Beta1_IRR[[2]])*10,na.rm=T)), col="darkblue",lwd=2,lty=3)
text(-140,1.00275,"E)",cex=1.2,xpd=TRUE)

plot(colMeans(HighConcurvity.Beta1_IRR[[1]], na.rm=T), type="n", ylim=c(0.990,1.010), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[1]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(HighConcurvity.Beta1_IRR[[3]])*10,na.rm=T)), col="lightskyblue1",lwd=2,lty=2)
lines(exp(colMeans(log(HighConcurvity.Beta1_IRR[[4]])*10,na.rm=T)), col="lightskyblue4",lwd=2,lty=3)
text(-140,1.015,"F)",cex=1.2,xpd=TRUE)


plot(colMeans(HighConcurvity.Beta2_IRR[[1]], na.rm=T), type="n", ylim=c(0.999,1.002), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[2]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(HighConcurvity.Beta2_IRR[[1]])*10,na.rm=T)), col="royalblue1",lwd=2,lty=2)
lines(exp(colMeans(log(HighConcurvity.Beta2_IRR[[2]])*10,na.rm=T)), col="darkblue",lwd=2,lty=3)
text(-140,1.00275,"G)",cex=1.2,xpd=TRUE)

plot(colMeans(HighConcurvity.Beta2_IRR[[1]], na.rm=T), type="n", ylim=c(0.990,1.010), ylab="Relative risk", xlab="Lag", bty="n",
     main="", xlim=c(0,750), xaxt="n")
axis(1,at=seq(0,800,by=100) )
abline(h=1)
lines(exp(Beta[[2]]*10),col="grey70", lwd=2, lty=1)
lines(exp(colMeans(log(HighConcurvity.Beta2_IRR[[3]])*10,na.rm=T)), col="lightskyblue1",lwd=2,lty=2)
lines(exp(colMeans(log(HighConcurvity.Beta2_IRR[[4]])*10,na.rm=T)), col="lightskyblue4",lwd=2,lty=3)
text(-140,1.015,"H)",cex=1.2,xpd=TRUE)

add_legend("bottom",cex=0.9,ncol=2,c("True","Not adjusted","NCS(doy,10df)+I(year)+\nNCS(week,5df)+NCS(month,5df)","NCS(t,4df/year)","NCS(t,10df/year)"),bty="n",
           col=c("grey","royalblue1","darkblue","lightskyblue1","lightskyblue4"),lwd=2,lty=c(1,2,3,2,3))
dev.off()




