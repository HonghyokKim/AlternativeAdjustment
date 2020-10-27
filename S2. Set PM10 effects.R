########True coefficients of PM10###########

T2length <- 731
T2 <- seq(T2length)

par(mfrow=c(3,4))
plot(a1<-1/((T2/4+exp(1))/1*log((T2/4+exp(1))/1))/6.433368*0.5877867*0.01)
plot(a2<-(0.22*(1-0.999)^(1)*0.999^T2)*0.01 )
plot(a3<-0.20*(1-0.99)^(1)*0.99^T2*0.01) 
plot(a4<-0.25*(1-0.9)^(1)*0.9^T2*0.01) 

#jpeg("%%%%FOLDERNAME%%%%/FIGURE_S1_Lagpattern.jpg",
#     width=15, height=15, pointsize=30,res=300,units="in")
layout(matrix(seq(8),ncol=2))
par(mar=c(4,4,2,1),mgp=c(2,1,0),oma=c(0,0,2,0))
plot(x=seq(0,730),y=a1<-1/((T2/4+exp(1))/1*log((T2/4+exp(1))/1))/6.433368*0.5877867*0.03/0.007249741*0.01,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Effect",cex=0.7)
text(x=-170, y=0.0016,"A",cex=1.2,xpd=T)
abline(h=0)
plot(x=seq(0,730),y=a2<- -(0.32*(1-0.92)^(1)*0.92^T2)*0.03/0.007249741*0.01,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Mortality displacement 1",cex=0.7)
text(x=-170, y=0.0003,"C",cex=1.2,xpd=T)
abline(h=0)
plot(x=seq(0,730),y=a3<- -0.1*(1-0.999)^(1)*0.999^T2*0.03/0.007249741*0.01,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Mortality displacement 2",cex=0.7)
text(x=-170, y=-0.0000012,"E",cex=1.2,xpd=T)
abline(h=0)
aa<-a1+a2+a3
plot(x=seq(0,730),y=aa<-(a1+a2+a3),type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Observable lag pattern 1",cex=0.7)
text(x=-170, y=0.000275,"G",cex=1.2,xpd=T)
abline(h=0)

plot(x=seq(0,730),y=b1<-1/((T2/4+exp(1))/1*log((T2/4+exp(1))/1))/6.433368*0.5877867*0.03/0.0007066529*0.001,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Effect",cex=0.7)
text(x=-170, y=0.0016,"B",cex=1.2,xpd=T)
abline(h=0)
plot(x=seq(0,730),y=b2<- -(0.25*(1-0.90)^(1)*0.90^T2)*0.03/0.0007066529*0.001,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Mortality displacement 1",cex=0.7)
text(x=-170, y=0.00029,"D",cex=1.2,xpd=T)
abline(h=0)
plot(x=seq(0,730),y=b3<- -0.35*(1-0.995)^(1)*0.995^T2*0.03/0.0007066529*0.001,type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Mortality displacement 2",cex=0.7)
text(x=-170, y=0.000024,"F",cex=1.2,xpd=T)
abline(h=0)
bb<-b1+b2+b3
plot(x=seq(0,730),y=bb<-(b1+b2+b3),type="l", bty="n", xlab="Lag",ylab="log-RR",xaxt="n")
axis(1,at=seq(0,730,by=100))
mtext(side=3,"Observable lag pattern 2",cex=0.7)
text(x=-170, y=0.000235,"H",cex=1.2,xpd=T)
abline(h=0)
dev.off()

Beta <- list("Type1"=matrix(NA,nrow=730+1,ncol=1),"Type2"=matrix(NA,nrow=730+1,ncol=1))

Beta[["Type1"]][,1]<-  aa
Beta[["Type2"]][,1]<-  bb
