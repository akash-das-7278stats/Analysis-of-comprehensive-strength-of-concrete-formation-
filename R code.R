data=read.csv("D:\\DU\\new project arijit\\slump_test.csv",header=TRUE) 
data=data[,-1]  
colnames(data)=c("cement","slag","fly","water","sp","ca","fa","slump","flow","cs") 
dt=data 
head(dt)
scaled=scale(dt[c("cement","slag","fly","water","sp","ca","fa")],method = "center")  
dt[c("cement","slag","fly","water","sp","ca","fa")]=scaled
head(dt)

pairs(dt[,-c(8,9)])
model1=lm(cs~.-(slump+flow),dt)
e=model1$residuals
b=coefficients(model1)
##PHASE !
par(mfrow=c(2,4))
for(i in 1:7)
{
  plot(dt[,i],dt$cs,xlab=colnames(data)[i],ylab="Compressive-strength")
}
#---Partial residual plot--
par(mfrow=c(2,4))
for(i in 1:7)
{
  plot(dt[,i],b[i+1]*dt[,i]+e,xlab=colnames(data)[i],ylab="Partial residual")
}
#-Added variable plot--
par(mfrow=c(2,4))
dt1=dt[-c(8,9)]

for(i in 1:7)
{
  plot(lm(dt1[,i]~.-(dt1$cs),data=dt1)$residuals,lm(cs~.-(dt1[,i]),data=dt1)$residuals,xlab=colnames(data)[i],ylab="y")
}
##Normality Test
par(mfrow=c(1,1))
qqnorm(e)
qqline(e)
library(dplyr)
shapiro.test(e)

##HETEROSCED
X=as.matrix(dt[,-(8:10)]) ; n=nrow(X)
X=cbind(rep(1,n),X)
H=X%*%solve(t(X)%*%X)%*%t(X)
b=c()
for(i in 1:n)
{
  bi=(e[i]^2)/(1-H[i,i])
}
yh=model1$fitted.values
plot(yh,b,xlab="fitted values",ylab="bi",main=" bi's vs fitted values")
plot(yh,e,xlab="fitted values",ylab="residuals",main="plot of residuals vs fitted values")
abline(h=0,col="red")

#bp test 
library(lmtest)
bptest(model1)

##godf q test
gqtest(model1)

#--Outlier detection--
A=dt[,-c(8,9)] 
r=c();t=c()
hat=H
p=ncol(X)
S=sqrt(sum(e^2)/(n-p))  
for(i in 1:n)
{ 
  r[i]=(e[i]/sqrt(1-hat[i,i]))/S 
  t[i]=r[i]*sqrt((n-p-1)/(n-p-r[i]^2))
}
plot(t)
abline(h=-2)
abline(h=2)
which(abs(t)>2)
#---covratio---
cr=c()
for(i in 1:n)
{ 
  cr[i]=1/(((((n-p-1)/(n-p))+((t[i]^2)/(n-p)))^p)*(1-hat[i,i]))
}
plot(cr,ylab="covration",main="Covratio")
abline(h=1.2,col="red")
abline(h=0.796,col="red")
#--cooks Distance--
D=c() 
for(i in 1:n) 
{ 
  D[i]=((r[i]^2)*hat[i,i])/(p*(1-hat[i,i])) 
} 
f=qf(0.10,p,n-p)
plot(D)
abline(h=f)
#--DFFITSi--
DS=c() 
for(i in 1:n)
{ 
  DS[i]=t[i]*sqrt(hat[i,i]/(1-hat[i,i]))
}
plot(DS)
abline(h=2*sqrt(p/n))
##---Testing hypothesis if i th obsn is an outlier-
which(abs(t)>qt(.975,n-p-1))
##phase 2
##truncating the outliers
outliers=c(8,14,22,49,60,87)
dt1=dt[-outliers,]
model2=lm(cs~.-(slump+flow),dt1); e=model2$residuals
##Normality Test
qqnorm(e)
library(dplyr)
shapiro.test(e)

#HETEROSCED
X1=as.matrix(dt1[,-(8:10)]) ;
n=nrow(X1)
X1=cbind(rep(1,n),X1)
H1=X1%*%solve(t(X1)%*%X1)%*%t(X1)
b=c() 
for(i in 1:n) 
{   b[i]=(e[i]^2)/(1-H1[i,i]) }
yh=model2$fitted.values 
plot(yh,b,xlab="fitted values",ylab="bi",main=" bi's vs fitted values")
library(lmtest)
bptest(model2)
gqtest(model2)

##AUTOCORR
e=model2$residuals ; n=length(e)
plot(e[-1],e[-n],cex=0.5,xlab="residuals",ylab="residuals",main="Plot of Residuals at lag 1")
abline(lm(e[-n]~e[-1]),col="red")
cor(e[-n],e[-1])
library(car)
durbinWatsonTest(model2,max.lag=1)
#---Outlier detection-
r1=c()
for(i in 1:n)
{
  r1[i]=(e[i]/sqrt(1-H1[i,i]))/(sum(e^2)/(n-p))
}
D1=c()
for(i in 1:n) 
{ 
  D1[i]=((r1[i]^2)*H1[i,i])/(p*(1-H[i,i])) 
} 
which(D1>f)
##MULTICOLLINEARITY
library(car)
vif(model2)
M=as.matrix(cor(X1[,-1]))
library(pracma)
cond(M)
#VIF of ca is highest 
model_3=lm(cs~.-(slump+flow+ca),dt1)
vif(model_3)
cond(as.matrix(cor(X1[,-c(1,7)])))# the condition became very small 
pairs(X1[,-1],pch=20)
round(cor(X1[,-1]),2)
### stepwise
intercept_only=lm(cs~1,data=dt1)
all_v= lm(cs~.-(slump+flow),data=dt1)
stepwise=step(intercept_only,direction="both",scope=formula(all_v)) ##put trace=0 to skip step details
summary(stepwise)
y_pred=stepwise$fitted.values
mad_step=mean(abs(y_pred-dt1$cs))
#lasso regression
library(caret)
X=dt1[,-(8:10)]
scaled=preProcess(X,method = "center") 
X=predict(scaled,X)
X=as.matrix(X)
dt1=as.data.frame(dt1)
library(pracma)
library(glmnet)
cs=dt1$cs
u1=cv.glmnet(X,cs,alpha = 1, family = 'gaussian') # alpha=1 means lasso
k1=0.02 ## how this value is selected(0.02912)?
lasso_reg=glmnet(X,cs,alpha = 1, family = 'gaussian', lambda = k1)
y_lasso=lasso_reg$a0+X%*%lasso_reg$beta
mad_lasso=mean(abs(y_lasso-dt1$cs))
m=0
for(i in 1:n)
{
 lasso_reg1=glmnet(X[-i,],cs[-i],alpha = 1, family = 'gaussian', lambda = k1)
 m=m+(cs[i]-lasso_reg1$a0-t(X[i,])%*%lasso_reg1$beta)^2
}
cv_lasso=m/n
cv_lasso
#----------- ridge ----------- ami korchi---------
u1=cv.glmnet(X,cs,alpha = 0, family = 'gaussian') # alpha=0 means ridge
k1=0.02 ## how this value is selected?
ridge_reg=glmnet(X,cs,alpha = 0, family = 'gaussian', lambda = k1)
y_ridge=ridge_reg$a0+X%*%ridge_reg$beta
mad_ridge=mean(abs(y_ridge-dt1$cs))
m=0
for(i in 1:n)
{
 ridge_reg1=glmnet(X[-i,],cs[-i],alpha = 0, family = 'gaussian', lambda = k1)
 m=m+(cs[i]-ridge_reg1$a0-t(X[i,])%*%ridge_reg1$beta)^2
}
cv_ridge=m/n
cv_ridge
#---------------------------------------------------
library(boot)
model_step1=glm(cs~.-(slag+slump+flow),data=dt)
cv=cv.glm(dt,model_step1)$delta
cv_step=cv[1]#cross validation error of model obtained by stepwise model selection

#comparison b/w the models

plot(cs,type="l",main="Original vs fitted response")
lines(y_pred,col="red")
lines(y_lasso,col="blue")
#--Least median square--
library(MASS)
dt=read.csv(file.choose(),header=TRUE)
dt=dt[,-1]
colnames(dt)=c("cement","slag","fly","water","sp","ca","fa","slump","flow","cs")
dt=as.matrix(dt[,-c(8,9)])
model_lms=lqs(dt[,8]~dt[,1]+dt[,2]+dt[,3]+dt[,4]+dt[,5]+dt[,6]+dt[,7],method = "lms")
y_lms=model_lms$fitted.values
plot(dt[,8],type="l",ylab="",main="Original vs fitted response")
lines(y_lms,col="green")
m1=0;coe=c()
for(i in 1:n)
{
  model_lms1=lqs(dt[-i,8]~dt[-i,1]+dt[-i,2]+dt[-i,3]+dt[-i,4]+dt[-i,5]+dt[-i,6]+dt[-i,7],method = "lms")
  coe=model_lms1$coefficients
  m1=m1+(cs[i]-coe[1]-t(X[i,])%*%lasso_reg$beta)^2
}
cv_lms=m1/n

