
library(mgcv)
library(splines)
library(tsModel)
library(dlnm)
library(lubridate)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

base_data <- read.csv("/Users/honghyokkim/Documents/Research/TS_overfitting/code/03312020/Code_Public/test_data_public.csv")
base_data$date <-as.Date(paste(base_data$day,base_data$month,base_data$year,sep='.'), format='%d.%m.%Y')

base_data$week_o <- week(base_data$date)
base_data$month_o <- base_data$month
for (ii in seq(2003,2013)) {
  base_data[base_data$year==ii,"week_o"] <- max(base_data[base_data$year==ii-1,"week_o"])+base_data[base_data$year==ii,"week_o"]
  base_data[base_data$year==ii,"month_o"] <- max(base_data[base_data$year==ii-1,"month_o"])+base_data[base_data$year==ii,"month_o"]
  
}
base_data$logPM10 <- log(base_data$PM10)

sim_dat <- base_data

tempknots <- quantile(sim_dat$ta_mean,c(10,50,90)/100,na.rm=TRUE)
tempbs <- crossbasis(sim_dat$ta_mean, lag=c(0,21),argvar=list(fun="bs",degree=2,knots=tempknots),arglag=list(knots=logknots(21, fun="bs", df=7)))

####Predict log_PM10
logPM10.fit <- (lm(logPM10 ~ ns(week_o,5)+ns(month_o,5)+ns(DAYOFYEAR,10)+as.factor(year)+as.factor(weekday_2)+as.factor(HOLIDAY), data=sim_dat))
logPM10.fit <- (lm(logPM10 ~ Lag(resid(logPM10.fit),k=1)+Lag(resid(logPM10.fit),k=2)+tempbs+ns(rh_mean,4)+
                     ns(week_o,5)+ns(month_o,5)+ns(DAYOFYEAR,10)+as.factor(year)+as.factor(weekday_2)+as.factor(HOLIDAY)+as.factor(inf_visit), data=sim_dat))
PM10.pred <- exp(predict(logPM10.fit))


sim_dat$week_o_1 <- ns(sim_dat$week_o,5)[,1]
sim_dat$week_o_2 <- ns(sim_dat$week_o,5)[,2]
sim_dat$week_o_3 <- ns(sim_dat$week_o,5)[,3]
sim_dat$week_o_4 <- ns(sim_dat$week_o,5)[,4]
sim_dat$week_o_5 <- ns(sim_dat$week_o,5)[,5]

sim_dat$month_o_1 <- ns(sim_dat$month_o,5)[,1]
sim_dat$month_o_2 <- ns(sim_dat$month_o,5)[,2]
sim_dat$month_o_3 <- ns(sim_dat$month_o,5)[,3]
sim_dat$month_o_4 <- ns(sim_dat$month_o,5)[,4]
sim_dat$month_o_5 <- ns(sim_dat$month_o,5)[,5]

####Predict Deaths
Y.fit <- glm(deaths~tempbs+ns(rh_mean,4)+week_o_1+week_o_2+week_o_3+week_o_4+week_o_5+
               month_o_1+month_o_2+month_o_3+month_o_4+month_o_5+
               ns(DAYOFYEAR,10)+as.factor(year)+
               as.factor(weekday_2)+as.factor(HOLIDAY)+as.factor(inf_visit),data=sim_dat, family=quasipoisson())
Y.pred <- exp(predict(Y.fit))
