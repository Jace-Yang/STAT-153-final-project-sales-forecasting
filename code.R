# 0 Data and packages loading
setwd("/Users/chichu/Desktop/Junior-1/Time series/Project")
rm(list=ls())
getwd()
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



# 1 Exploratory Data Analysis

## Data Transform
sales_raw <- read.csv("sales.csv")
date = (mdy(sales_raw[,2]))
sales_value = sales_raw[,3]
sales = zooreg(sales_value, order.by=date, start = date[1])
sales_ts = ts(sales, frequency = 7)
sales_frame = data.frame(time = date, sales = sales)

# 1.1 Original plot
ggthemr_reset()
theme_set(theme_gray(base_size = 16))

plot_original = 
  ggplot() + 
    aes(x = date, y = sales_value, colour = sales_value)+
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
plot_season = ggseasonplot(ts(daily2monthly(sales,FUN=mean,na.rm=TRUE),frequency = 12,start=c(2015,1)))+ylab("Monthly Average")+
    ggtitle("Seasonal Plot")+theme(plot.title = element_text(
      colour = "black",size = 18, 
      face = "bold", hjust = 0.5))

## 1.3 Box-plot
sales_box <- daily2monthly(sales,FUN=mean) %>% ts(frequency = 12)
Month <- factor(x= as.vector(cycle(sales_box)), levels = 1:12, labels = c("Jan", "Feb", "Mar", "Apr", 
                                                                      "May", "Jun", "Jul", "Aug",
                                                                      "Sep", "Oct", "Nov", "Dec"))
plot_box = ggplot(data.frame(Sales = as.vector(sales_box), Month))+ 
  geom_boxplot(aes(x = Month, y = Sales, group = Month))+ylab("Monthly Average Sales")+
  ggtitle("Box-Plot")+ theme(plot.title = element_text(size = 18, face="bold", hjust = 0.5))

## Summary
ggarrange(plot_original,ggarrange(plot_season,plot_box,ncol=2,nrow=1),ncol=1,nrow=2)


## CV of Models

### Function of plotting
plot_predict = function(model_summary, name){
  theme_set(theme_gray(base_size = 16))
  ggplot(model3_summary)+
    geom_line(aes(x=date,y=sales_value, colour = "Actual"), size=0.8,na.rm=TRUE) +
    geom_line(aes(x=date,y=Fitted, colour = "Fitted"), size=0.5,na.rm=TRUE)+
    geom_line(aes(x=date,y=Predict, colour = "Predict"), size=0.5,na.rm=TRUE)+
    scale_color_manual(values=c("#e87d71", "#56B4E9", "#1f78b4"))+
    ggtitle(paste("In-sample and Out-sample fitting VS real data in" ,name))+
    labs(colour = "Data type")+ylab("Sales")+xlab("Years")+
    theme(plot.title = element_text(face="bold", hjust = 0.5,vjust=-1.9,size = 19,margin=margin(0,0,30,0)))
}
### CV of Model 1

### CV of Model 2

### CV of Model 3

N=1200-365
a=XSIN_Trend(sales_value[1:N],pred_num=length(sales_value)-N,WarnOnly = T)


model3_summary = data.frame(date=date,
       Actual = sales_value,
       Fitted = c(a[1:N],rep(NA,length(sales_value)-N)),
       Predict = c(rep(NA, N),a[(N+1):length(sales_value)]))

### CV summary
plot_model1_summary = plot_predict(model1_summary,"Model 3")
plot_model2_summary = plot_predict(model2_summary,"Model 3")
plot_model3_summary = plot_predict(model3_summary,"Model 3")
ggarrange(plot_model1_summary,plot_model1_summary,plot_model1_summary,ncol=1,nrow=3)




## ggsarima
ggthemr_reset()
ggsarima = function(residuals,name){
  residuals = as.vector(residuals)
  plot_residuals = ggplot()+
    aes(1:length(residuals), residuals) +
    geom_line(color="#1f78b4")+ xlab("Index of Residuals")+
    labs(title = paste("Residuals of" , name))+
    theme(plot.title = element_text(size = 19, hjust = 0.5, face="bold"),
          axis.title=element_text(size=13))
    
  
  plot_ACF = ggAcf(residuals,lag.max = 365)+ggtitle(paste("ACF of" ,name,"'s Residuals"))+
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))

  plot_PACF = ggAcf(residuals,lag.max = 365)+ggtitle(paste("PACF of" ,name,"'s Residuals"))+ 
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))

  plot_gg = ggqqplot(as.vector(residuals),ggtheme = theme_gray())+
    ggtitle(paste("Q-Q Plot of" ,name,"Std Residuals"))+ 
    theme(plot.title = element_text(face="bold", hjust = 0.5),
          axis.title=element_text(size=13))
  
  ljungbox <- function(i, resi=residuals) {
    Box.test(resi, i, type = "Ljung-Box")$p.value
  }
  lag = 1:20
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
model =  auto.arima(sales_value)
residuals = model$residuals
name = "fuck"
ggsarima(residuals,"Model 1")



ggsarima_simple = function(residuals,name){
  residuals = as.vector(residuals)
  plot_residuals = ggplot()+
    aes(1:length(residuals), residuals) +
    geom_line(color="#1f78b4")+ xlab("Index of Residuals")+
    labs(title = paste("Residuals of" , name))+
    theme(plot.title = element_text(size = 19, hjust = 0.5, face="bold"),
          axis.title=element_text(size=13))
  
  plot_ACF = ggAcf(residuals,lag.max = 365)+ggtitle(paste("ACF of" ,name,"'s Residuals"))+
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))
  
  plot_PACF = ggAcf(residuals,lag.max = 365)+ggtitle(paste("PACF of" ,name,"'s Residuals"))+ 
    theme(plot.title = element_text(face="bold", hjust = 0.5,size = 15),
          axis.title=element_text(size=13))
  
  ggarrange(plot_residuals,
            ggarrange(plot_ACF,plot_PACF,ncol=2,nrow=1),
            ncol=1,nrow=2)
}

model =  auto.arima(sales_value)
residuals = model$residuals
name = "fuck"
ggsarima(residuals,"Model 1")



# What does tsdiag look like?
fit <- arima(lh, c(1,0,0))
tsdiag(fit)

# Need to extract the calculations from the existing
# plotting code.  I started by looking at the contents of
# stats:::tsdiag.Arima

# Standardised residuals


# To figure out the confidence bands, we need to dig into the code
# for stats:::plot.acf.  Unfortunately this isn't available in its own
# method, so I'm not entirely sure I've extracted it correctly
ci <- 0.95
clim <- qnorm((1 + ci)/2) / sqrt(acf$n.used)

qplot(lag, acf, data = acf_df, yend = 0, xend = lag, geom="segment") +
  geom_hline(colour = "grey50") +
  geom_point() +
  geom_hline(yintercept = c(-clim, clim), colour = "darkblue")


# Finally the code for the Ljung-Box statistic


size = 15


data("ToothGrowth")
df <- ToothGrowth
df$dose <- as.factor(df$dose)
# Create some plots
# ::::::::::::::::::::::::::::::::::::::::::::::::::
# Box plot
bxp <- ggboxplot(df, x = "dose", y = "len",
                 color = "dose", palette = "jco")


library(tidyr)
gather

mtcars_longer = mtcars %>% gather(attribute, value, cyl:carb) 
