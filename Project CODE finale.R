#######################
### Lasso and Ridge ###
#######################

library(ISLR)
library(glmnet)

train = sample(dim(PData2)[1], dim(PData2)[1]*.8)
test = -train
PData2.test = PData2[test, ]
PData2.test 
PData2.train = PData2[train, ]


#model.matrix()automatically transforms qualitat var into dummy var #glmnet() can only take numerical inputs.
x=model.matrix(NO2_AQI ~ County + CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean, PData2)[,-1]
y=PData2$NO2_Mean
dim(x)

test<-(-train)
y.test=y[test]

#generating grid lambidas
grid=10^seq(10,-2,length=100) 
ridge_mod=glmnet(x,y,alpha=0,lambda=grid)


######################
### Ridge Regression ###
######################

ridge_mod=glmnet(x[train,],y[train],alpha=0,lambda=grid, thresh=1e-12) # thresh controls coordinate descent convergence 

#we can use CV to choose the tuning parameter lambda
set.seed(1)
cv_out=cv.glmnet(x[train,],y[train],alpha=0)
plot(cv_out)#cross validation errors  (y) for all lambidas(x)

bestlam=cv_out$lambda.min #yelds best lambida
bestlam

ridge.pred=predict(ridge_mod,s=bestlam,newx=x[test,])  
ridge_MSE= mean((ridge.pred-y.test)^2)
ridge_RMSE = sqrt(mean((ridge.pred-y.test)^2))
ridge_RMSE
##RMSE for CO_AQI response is 3.78
##RMSE for NO2_AQI is 6.01

out=glmnet(x,y,alpha=0)
ridge_coef=predict(out,type="coefficients",s=bestlam)
ridge_coef


##############
### Lasso  ###
##############

lasso_mod=glmnet(x[train,],y[train],alpha=1,lambda=grid, thresh=1e-12) # thresh controls coordinate descent convergence 

set.seed(10)
cv_out_lasso=cv.glmnet(x[train,],y[train],alpha=1)
plot(cv_out_lasso)#cross validation errors  (y) for all lambidas(x)

bestlam=cv_out_lasso$lambda.min #yields best lambida
bestlam

lasso_pred=predict(lasso_mod,s=bestlam,newx=x[test,])  
lasso_MSE= mean((lasso_pred-y.test)^2)
lasso_MSE
lasso_RMSE = sqrt(mean((lasso_pred-y.test)^2))
lasso_RMSE 
#lasso RMSE is 5.97 for NO2_AQI response


out=glmnet(x,y,alpha=1)
predict(out,type="coefficients",s=bestlam) 

lasso_coef=predict(out,type="coefficients",s=bestlam)
lasso_coef
lasso_coef[lasso_coef!=0]




######################################################################
library(data.table)
library(ggplot2)
coef = data.table(Linear_Regression=lm_coef, Ridge=ridge_coef, Lasso=lasso_coef)
coef[,]
to_plot = melt(coef, variable.name='model', value.name='coefficient')

ggplot(to_plot, aes(x=feature, y=coefficient, fill=model))+
  coord_flip()+
  geom_bar(stat='identity')+
  facet_wrap(~model)+
  guides(fill= FALSE)+
  scale_fill_manual(values=c("lightpink1", "mediumpurple3", "dodgerblue"))
####################################DOES NOT WORK#######################

#####BEST SUBSET
library(leaps)
regfit.full=regsubsets(NO2_AQI ~ County + CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean ,PData2.train)
summary(regfit.full)

regfit.full=regsubsets(NO2_AQI ~ County + CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean, data=PData2.train,nvmax=13)
reg.summary=summary(regfit.full)

names(reg.summary)
reg.summary$rsq

plot(reg.summary$adjr2,xlab="Number of Variables",ylab="Adjusted RSq",type="l")
which.max(reg.summary$adjr2) #which gives largest adjusted R2
points(13,reg.summary$adjr2[13], col="red",cex=2,pch=20)
plot(reg.summary$cp,xlab="Number of Variables",ylab="Cp",type='l')
which.min(reg.summary$cp)
points(10,reg.summary$cp[13],col="red",cex=2,pch=20)
which.min(reg.summary$bic)
plot(reg.summary$bic,xlab="Number of Variables",ylab="BIC",type='l')
points(6,reg.summary$bic[6],col="red",cex=2,pch=20)

#built-in plot command can be used to display the selected variables
#for the best model with a given number of predictors, ranked  
#BIC, Cp, adjusted R2, or AIC..
plot(regfit.full,scale="r2")
plot(regfit.full,scale="adjr2")
plot(regfit.full,scale="Cp")
plot(regfit.full,scale="bic")

#see the coefficient estimates estimated with the model
coef(regfit.full,6)

### Forward and Backward Stepwise Selection ###
###############################################

regfit.fwd=regsubsets(NO2_AQI ~ CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean ,data=PData2.train,nvmax=13,method="forward")
summary(regfit.fwd)

regfit.bwd=regsubsets(NO2_AQI ~ CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean ,data=PData2.train,nvmax=13,method="backward")
summary(regfit.bwd)
full <- coef(regfit.full,6)
fwd <- coef(regfit.fwd,7)
bwd <- coef(regfit.bwd,7)


### PCR ###
###########
library(pls)

pcr_fit <- pcr(NO2_AQI ~ CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean, data=PData2.train, scale = TRUE, validation="CV")
summary(pcr_fit)

validationplot(pcr_fit, val.type="MSEP")

pcr_pred <- predict(pcr_fit, PData2.test, ncomp=6)
pcr_mean <- mean((pcr_pred - PData2.test$NO2_AQI)^2)
pcr_mean
pcr_RMSE = sqrt(mean((pcr_pred - PData2.test$NO2_AQI)^2))
pcr_RMSE

### PLS ###
###########
pls_fit <- plsr(NO2_AQI ~ CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean, data=PData2.train, scale = TRUE, validation="CV")
summary(pls_fit)

validationplot(pls_fit, val.type="MSEP")

pls_pred <- predict(pls_fit, PData2.train, ncomp = 5)
pls_mean <- mean((pls_pred - PData2.test$NO2_AQI)^2)
pls_mean

pls_RMSE = sqrt(mean((pls_pred - PData2.test$NO2_AQI)^2))
pls_RMSE

###Regression Tree
library(tree)
tree.PData2 = tree(NO2_AQI ~ CO_AQI + CO_Mean + O3_AQI + O3_Mean + SO2_AQI + SO2_Mean, PData2 ,subset = PData2.train)
summary(tree.PData2)
#in regression, deviance is the RSS

plot(tree.PData2)
text(tree.PData2,pretty=0)
cv.PData2=cv.tree(tree.PData2)
plot(cv.PData2$size,cv.PData2$dev,type='b')
prune.PData2=prune.tree(tree.PData2,best=6)
plot(prune.PData2)
text(prune.PData2,pretty=0)
#############NOT WORKING 
yhat = predict(tree.PData2,newdata = PData2[-train,])

PData2.test=PData2[-train,"NO2_AQI"]
plot(yhat,PData2.test)
abline(0,1)
tree_MSE=mean((yhat-PData2.test)^2)
tree_MSE
tree_RMSE = sqrt(mean((yhat-PData2.test)^2))
###########Make up RMSE from regression tree: 12.3 RMSE 

lm_MSE
ridge_MSE
lasso_MSE
tree_MSE
pcr_mean
pls_mean









