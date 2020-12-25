# setwd("~/Rs/TS_Project")
setwd("/Users/yaochen/Desktop/UCB_Study/Time_Series/Project")
rm(list=ls())
dev.off()

###################### load package and data #############################
library(lubridate)
library(astsa)
library(forecast)
library(dplyr)
library(lubridate)
library(forecast)
library(zoo)
library(xts)
library(ggplot2)
library(ggthemes)
library(ggthemr)
library(ggThemeAssist)
library(ggfortify)
library(hydroTSM)
library(ggpubr)
library(ggDiagnose)
sales_raw <- read.csv("sales.csv")
sales_value <- sales_raw[,3]
start_date=as.Date("2015-01-01")

####################### Exploratory Analysis #############################
## Data Transform
time = (mdy(sales_raw[,2]))
sales_value = sales_raw[,3]
sales = zooreg(sales_value, order.by=time, start = time[1])
sales_ts = ts(sales, frequency = 7)
sales_frame = data.frame(time = time, sales = sales)

# 1.1 Original plot
ggthemr_reset()
theme_set(theme_gray(base_size = 16))

plot_original = 
  ggplot() + 
  aes(x = time, y = sales_value, colour = sales_value)+
  geom_line(size = 0.4) +
  labs(colour = "Sales")+
  ggtitle("Original Sales")  + xlab("Year")+ylab("Sales")+
  theme(plot.title = element_text(
    colour = "black",size = 22, 
    face = "bold", hjust = 0.5,vjust = -0.5,
    margin=margin(0,0,15,0)),
    legend.title = element_text(size = 13))
scale_color_gradient(low = "#319cd0", high = "#0E294B") 

## 1.2 Seasonal Plot
plot_season = ggseasonplot(ts(daily2monthly(sales,FUN=mean,na.rm=TRUE),
                              frequency = 12,start=c(2015,1)))+
  ylab("Monthly Average")+
  ggtitle("Seasonal Plot")+theme(plot.title = element_text(
    colour = "black",size = 18, 
    face = "bold", hjust = 0.5))

## 1.3 Box-plot
sales_box <- daily2monthly(sales,FUN=mean) %>% ts(frequency = 12)
Month <- factor(x= as.vector(cycle(sales_box)), levels = 1:12, 
                labels = c("Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", 
                           "Aug","Sep", "Oct", "Nov", "Dec"))
plot_box = ggplot(data.frame(Sales = as.vector(sales_box), Month))+ 
  geom_boxplot(aes(x = Month, y = Sales, group = Month))+
  ylab("Monthly Average Sales")+
  ggtitle("Box-Plot")+ 
  theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5))

## Summary
ggarrange(plot_original,ggarrange(plot_season,plot_box,ncol=2,nrow=1),
          ncol=1,nrow=2)


########################### Model1 Linear ################################
# Fit the deterministic part
## create the data frame 
sales0 <- sales_value
date <- start_date+1:length(sales0)-1
time <- 1:length(sales0)
salesDf <- data.frame(time, date, sales0)

## delete all the data before 2016-03-01
salesDfNoLeap <- salesDf[-(1:425),]

## extract the time series sales
sales_ts <- ts(salesDfNoLeap[,"sales0"])

## variance stablization transformation
lambda <- BoxCox.lambda(sales_ts+10)
trans_sales <- BoxCox(sales_ts+10, lambda = lambda)

## create indicator variables for 4th July, high seasons
month <- salesDfNoLeap[,"date"] %>% month()
weekday <- salesDfNoLeap[,"date"] %>% wday() %>% as.factor()
day <- salesDfNoLeap[,"date"] %>% day()
Indicator4th <- day == 4 & month == 7
IndicatorMonth <- month %in% c(5,6,7,8,9)

## form the new data frame
salesDfNoLeap <- cbind(salesDfNoLeap,trans_sales = as.numeric(trans_sales), 
                       weekday, month, Indicator4th, IndicatorMonth)
salesDfNoLeap$time <- 1:nrow(salesDfNoLeap)

## fit the linear model
model1 <- lm(trans_sales~Indicator4th + IndicatorMonth + weekday + 
               time + time*IndicatorMonth, 
             data = salesDfNoLeap)

# Fit the residuals
s1 <- sarima(model1$residuals, 7,0,0,0,0,0)

ggsarima(s1$fit$residuals,"Model 1: AR(7)")


######################### Model2 Difference ##############################
# Difference
diffTsNoleap <- diff(diff(trans_sales, lag = 7, differences = 2), lag = 365)

# Check the data after differencing
plot.ts(diffTsNoleap)
acf2(diffTsNoleap)

# Fit the residuals
s2 = sarima(diffTsNoleap, 0,0,7,0,0,1,7)
ggsarima(s2$fit$residuals,"Model 2: ARIMA(0,0,7)(0,0,1)[7]")


######################### Model3 Parametric ##############################
# Fit the deterministic part

ND=length(sales_value)
t_in  = 1:ND
start_date=as.Date("2015-01-01")
library(lubridate)
library(astsa)
library(forecast)

sdate=start_date+t_in-1
if(ND<365) stop("Less than 1 year,periodical trend cannot be fitted")

sales=sales_value/((t_in^0.33) + 10) 

July4=(day(sdate)==4&month(sdate)==7)
weekd=as.factor(wday(sdate))
model0 = lm( sales~ weekd+July4:t_in) 

plain=c(1:110,290:460,660:800,1030:1180,1400:1530)[1:ND]
ptrend=lm(model0$residuals[plain]~t_in[plain]) 

deTreDum=model0$residuals-ptrend$coefficients[1]-ptrend$coefficients[2]*t_in

s=pmax(0,sin(t_in/365*2*pi-pi*7/12))
trend1=lm(deTreDum~s:t_in+s)
trend5=nls(deTreDum~pmax(0,(a*t_in^2+b)*sin(t_in/365*2*pi-phi)+c*t_in+d),
           start = list(a=trend1$coefficients["s:t_in"], 
                        b=trend1$coefficients["s"],
                        c=0,d=0,phi=pi*7/12),
           control=nls.control(warnOnly =T)) 

YDa=deTreDum-predict(trend5)

# plot trend and residuals for model3
plot_trend = ggplot() + 
  aes(x = date, y = deTreDum, colour = sales_value)+
  geom_line(size = 0.4) +
  geom_line(aes(x=date,y=predict(trend5)),color="red",size=1.15)+
  labs(colour = "Sales")+
  ggtitle("Fitted trend for pre-detrend data")  + xlab("Year")+ylab("Sales")+
  theme(plot.title = element_text(
    colour = "black",size = 15, 
    face = "bold", hjust = 0.5,vjust = -0.5),
    axis.title=element_text(size=15),
    legend.title = element_text(size = 13))+
    scale_color_gradient(low = "#319cd0", high = "#0E294B")

plot_residuals_trend = ggplot()+
  aes(1:length(YDa), YDa) +
  geom_line(color="#1f78b4")+ xlab("Index of Residuals")+ylab("Residuals")+
  labs(title = "Residuals of finale de-trend Data" )+
  theme(plot.title = element_text(size = 15, hjust = 0.5,vjust = -0.5 ,
                                  face="bold"),
        axis.title=element_text(size=15))
ggarrange(plot_residuals_trend,plot_trend,nrow=1,ncol=2)+
  ggtitle("Model 3 Trend and Residuals")+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face="bold"))

# Fit the residuals
s3 = sarima(YDa,5,0,0,0,1,1,7) 

# ggsarima for residuals
ggsarima = function(residuals,name){
  residuals = as.vector(residuals)
  plot_residuals = ggplot()+
    aes(1:length(residuals), residuals) +
    geom_line(color="#1f78b4")+ xlab("Index of Residuals")+
    labs(title = paste("Residuals of" , name))+
    theme(plot.title = element_text(size = 19, hjust = 0.5, face="bold"),
          axis.title=element_text(size=13))
  
  
  plot_ACF = ggAcf(residuals,lag.max = 365)+
    ggtitle(paste("ACF of" ,name,"'s Residuals"))+
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))
  
  plot_PACF = ggAcf(residuals,lag.max = 365)+
    ggtitle(paste("PACF of" ,name,"'s Residuals"))+ 
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))
  
  plot_gg = ggqqplot(as.vector(residuals),ggtheme = theme_gray())+
    ggtitle(paste("Q-Q Plot of" ,name,"Std Residuals"))+ 
    theme(plot.title = element_text(face="bold", hjust = 0.5),
          axis.title=element_text(size=13))
  
  ljungbox <- function(i, resi=residuals) {
    Box.test(resi, i, type = "Ljung-Box")$p.value
  }
  lag = 1:50
  Pvalue <- sapply(lag,ljungbox)
  plot_Ljungbox = 
    ggplot() + 
    geom_point(aes(x=lag, y=Pvalue),color="black",pch=1)+xlab("Lag")+
    geom_hline(yintercept = 0.05, colour = "#3962f3",linetype=2)+ylim(0,1)+
    theme(plot.title = element_text(size = 15, hjust = 0.5, face="bold"),
          axis.title=element_text(size=13))+ 
    labs(title = paste("P-values for Ljung-Box in", name))
  
  ggarrange(plot_residuals,
            ggarrange(plot_ACF,plot_PACF,ncol=2,nrow=1),
            ggarrange(plot_gg,plot_Ljungbox,ncol=2,nrow=1),
            ncol=1,nrow=3)
}

residuals = s3$fit$residuals
ggsarima(residuals,"Model 3: ARIMA(5,0,0)(0,1,1)[7]")

######################### Cross Validation ###############################
CrossVali <- function(Tdata,Nfun){
  Start_t=365+366
  Pre=1
  MSE=foreach(i=Start_t:length(Tdata),.combine='c') %dopar% {
    mean((Tdata[i-1+1:Pre]-tail(Nfun(sales0=Tdata[1:(i-1)],
                                     pred_num=Pre),N=Pre))^2)
  }
  return(sum(MSE))
}

# model1 forecast function(sales is a vector of raw sales data)
model1For <- function(sales0=sales_value, pred_num = 1){
  start_date=as.Date("2015-01-01")

  # create the data frame 
  date <- start_date+1:length(sales0)-1
  time <- 1:length(sales0)
  salesDf <- data.frame(time, date, sales0)
  
  # delete all the data before 2016-03-01
  salesDfNoLeap <- salesDf[-(1:425),]
  
  # extract the time series sales
  sales_ts <- ts(salesDfNoLeap[,"sales0"])
  
  # variance stablization transformation
  lambda <- BoxCox.lambda(sales_ts+10)
  trans_sales <- BoxCox(sales_ts+10, lambda = lambda)
  
  # create indicator variables for 4th July, high seasons
  month <- salesDfNoLeap[,"date"] %>% month()
  weekday <- salesDfNoLeap[,"date"] %>% wday() %>% as.factor()
  day <- salesDfNoLeap[,"date"] %>% day()
  Indicator4th <- day == 4 & month == 7
  IndicatorMonth <- month %in% c(5,6,7,8,9)
  
  # form the new data frame
  salesDfNoLeap <- cbind(salesDfNoLeap,trans_sales = as.numeric(trans_sales), 
                         weekday, month, Indicator4th, IndicatorMonth)
  salesDfNoLeap$time <- 1:nrow(salesDfNoLeap)
  
  # fit the linear model
  model1 <- lm(trans_sales~Indicator4th + IndicatorMonth + weekday + 
                 time + time*IndicatorMonth, 
               data = salesDfNoLeap)
  
  # fit the residuals
  resModel1 <- sarima(model1$residuals, 7,0,0,0,0,0)
  
  # forecast the residuals
  resFor <- sarima.for(model1$residuals, n.ahead = pred_num, 7,0,0,0,0,0)$pred
  
  # forecast the linear model
  preDate <- tail(salesDfNoLeap$date,1)+1:pred_num
  newDate <- c(salesDfNoLeap$date, preDate)
  newMonth <- month(newDate)
  newWeek <- wday(newDate) %>% as.factor()
  newDay <- day(newDate)
  newIndicator4th <- newDay == 4 & newMonth == 7
  newIndicatorMonth <- newMonth %in% c(5,6,7,8,9)
  newTime <- 1:length(newDate)
  predictDF <- data.frame(Indicator4th = newIndicatorMonth, 
                          IndicatorMonth = newIndicatorMonth,
                          weekday = newWeek, 
                          time = newTime)
  predTran <- predict.lm(model1, newdata = predictDF) + 
    c(model1$residuals, resFor)
  
  # transform the predict value into raw data
  predRaw <- (predTran*lambda +1)^(1/lambda)
  return(predRaw)
}

# cross validation for model2
forecast_residual <- c()
forecast_value <- c()
train_set <- c()
error <- 0
for (i in 381:1217) {
  train_set <- diffTsNoleap[1:(i-379)]
  forecast_residual[i] <- sarima.for(train_set,n.ahead = 1,0,0,7,0,0,1,7)$pred
  forecast_value[i] <- forecast_residual[i] + 2*trans_sales[i-7] - trans_sales[i-14]
  + trans_sales[i-365] - 2*trans_sales[i-372] + trans_sales[i-379]
  error <- error + ((forecast_value[i] - sales_ts[i])^2)[1]
}
error

# model3 forecast function(sales is a vector of raw sales data)
XSIN_Trend<-function(sales0=sales_value,pred_num=1,Peak=2,
                     tin_times=0.33,tin_int=10,WarnOnly=T)
{
  ND=length(sales0)
  t_in  = 1:ND
  start_date=as.Date("2015-01-01")
  library(lubridate)
  library(astsa)
  library(forecast)
  
  sdate=start_date+t_in-1
  if(ND<365) stop("Less than 1 year,periodical trend cannot be fitted")
  
  sales=sales0/((t_in^tin_times) + tin_int) ## 1.
  
  July4=(day(sdate)==4&month(sdate)==7)
  weekd=as.factor(wday(sdate))
  model0 = lm( sales~ weekd+July4:t_in) ## 2.
  
  plain=c(1:110,290:460,660:800,1030:1180,1400:1530)[1:ND]
  ptrend=lm(model0$residuals[plain]~t_in[plain]) ## 3.
  deTreDum=model0$residuals-ptrend$coefficients[1]-ptrend$coefficients[2]*t_in
  
  s=pmax(0,sin(t_in/365*2*pi-pi*7/12))
  trend1=lm(deTreDum~s:t_in+s)
  trend5=nls(deTreDum~pmax(0,(a*t_in^Peak+b)*sin(t_in/365*2*pi-phi)+c*t_in+d),
             start = list(a=trend1$coefficients["s:t_in"], 
                          b=trend1$coefficients["s"],
                          c=0,d=0,phi=pi*7/12),
             control=nls.control(warnOnly =WarnOnly)) ## 4.
  YDa=deTreDum-predict(trend5)
  
  # s1 = sarima(Yd_,0,0,7,0,0,1,7)
  # acf2(s1$fit$residuals,max.lag = 800)
  # 
  # s1 = sarima(Yd_,0,0,0,0,0,1,7)
  # acf2(s1$fit$residuals,max.lag = 800)
  # 
  # s1 = sarima(Yd_,5,0,0,0,0,1,7)
  # acf2(s1$fit$residuals,max.lag = 800)
  s1 = sarima(YDa,5,0,0,0,1,1,7) ## 5.
  
  
  if(pred_num<0) stop("pred_num>=0!!!")
  t_in=1:(ND+pred_num)
  sdate=start_date+t_in-1
  July4=(day(sdate)==4&month(sdate)==7)
  weekd=as.factor(wday(sdate))
  
  if(pred_num!=0)  Fitted=c(YDa-s1$fit$residuals,
                      as.numeric(sarima.for(YDa,pred_num,5,0,0,0,1,1,7)$pred))
  else Fitted=YDa-s1$fit$residuals
  Fitted=Fitted+predict(trend5,newdata = list(t_in=1:(ND+pred_num)))
  Fitted=Fitted+ptrend$coefficients[1]+ptrend$coefficients[2]*t_in
  Fitted=Fitted+predict(model0,newdata = list(weekd=weekd,July4=July4,t_in=t_in))
  Fitted[July4]=Fitted[July4]-mean(Fitted) # Adjusted lm
  Fitted=Fitted*((t_in^tin_times) + tin_int)
  
  # dev.off()
  return(Fitted)
}

# cross validation for model 1 and 3
a=tail(XSIN_Trend(pred_num=10,WarnOnly = T),10)
plot(1543:1642,tail(c(sales_value),100),type="l",xlim = c(1543,1652),ylim = c(0,260))
lines(1643:1652,a,col="red")

a=tail(model1For(pred_num=10),10)
plot(c(sales_value),type="l",xlim = c(0,length(sales_value+10)),ylim=c(0,300))
lines(1643:1652,a,col="red")

jpeg("TRASH.jpg")
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
TIME<-system.time({
  SSE=CrossVali(sales_value,XSIN_Trend)
})
stopCluster(cl)
dev.off()

jpeg("TRASH.jpg")
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
TIME<-system.time({
  SSE=CrossVali(sales_value,model1For)
})
stopCluster(cl)
dev.off()


# Plots of Cross Validation(red line is in-sample prediction, blue line is out-sample)
a=XSIN_Trend(pred_num=10,WarnOnly = T)
plot(c(sales_value),type="l",xlim = c(0,length(sales_value+10)))
lines(a,col="red")


## CV of Models
model_summary$date

### Function of plotting
plot_predict = function(model_summary, name){
  theme_set(theme_gray(base_size = 16))
  ggplot(model_summary)+
    geom_line(aes(x=date,y=Actual, colour = "Actual"), size=0.8,na.rm=TRUE) +
    geom_line(aes(x=date,y=Fitted, colour = "Fitted"), size=0.5,na.rm=TRUE)+
    geom_line(aes(x=date,y=Predict, colour = "Predict"), size=0.5,na.rm=TRUE)+
    scale_color_manual(values=c("#e87d71", "#56B4E9", "#1f78b4"))+
    ggtitle(paste("In-sample and Out-sample fitting VS real data in" ,name))+
    labs(colour = "Data type")+ylab("Sales")+xlab("Years")+
    theme(plot.title = element_text(face="bold", hjust = 0.5,vjust=-1.9,
                                    size = 19,margin=margin(0,0,30,0)))
}

### CV of Model 3

N=1200
a=XSIN_Trend(sales_value[1:N],pred_num=length(sales_value)-N,WarnOnly = T)

model3_summary = data.frame(date=date,
                            Actual = sales_value,
                            Fitted = c(a[1:N],rep(NA,length(sales_value)-N)),
                            Predict = c(rep(NA, N),a[(N+1):length(sales_value)]))

### CV summary
plot_model3_summary = plot_predict(model3_summary,"Model 3")


