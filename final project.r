install.packages('data.table')
install.packages('MASS')
install.packages('VGAM')
install.packages("ggplot2")
install.packages("ggrepel")
install.packages('ggpubr')
install.packages('ggExtra')
install.packages('stargazer')
install.packages('corrplot')
install.packages('olsrr')
install.packages('eply')

library(data.table)
library(MASS)
library(car)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggExtra)
library(stargazer)
library(nlme)
library(corrplot)
library(Matrix)
library(olsrr)
library(eply)
setwd("C:/Users/Han/Desktop/Box Sync/Stat 501/Final Project")

dt_1 = fread('final project dataset 1.txt')
dt_2 = fread('final project dataset 2.txt')

dt = merge(dt_1,dt_2,by = 'State')

#[1] 0
## No NA value

#### 1.	Find an appropriate linear regression model with M as the dependent variables.
#### Among the independent variables include MA, D, PL, S, B, HT, UR, CR, HS, INC, PL, VT, and UR. 
#### Add and subtract variables as you think fit.
#### Notice that you may need to detect influential points, multicollineaity, also do necessary transformation if needed.

dt1 = dt[,c('State','M','MA','D','PL','S','B','HT','UR','CR','HS','INC','VT')]
dt1[,2:13] = dt1[,2:13][,lapply(.SD,as.numeric)]
# standardize all X
dt1[,3:13] = dt1[,3:13][,lapply(.SD,scale)]

# orignal EDA
temp <- ggplot(dt1, aes(x=MA, y=M)) + geom_point() + theme_classic() + geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x11 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=D, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x12 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=PL, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x13 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=S, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x14 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=B, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x15 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=HT, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x16 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=UR, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x17 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=CR, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x18 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=HS, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x19 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=INC, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x110 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1, aes(x=VT, y=M)) + geom_point() + theme_classic()+ geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
x111 = ggMarginal(temp,type = 'histogram', fill = "cyan")
temp <- ggplot(dt1[-33,], aes(x=MA, y=M)) + geom_point() + theme_classic() + geom_smooth(method = 'lm')+ theme(text = element_text(size=18, face = 'bold'))+
        labs(x = 'influential pt temporarily removed MA')
x112 = ggMarginal(temp,type = 'histogram', fill = "cyan")
ggarrange(x11,x112,x12,x13,x14,x15,x16,x17,x18,x19,x110,x111,ncol = 4,nrow = 3)

tr <- function (m)
{
   total_sum <- 0
   if(is.matrix(m))
   {
      row_count <- nrow(m)
      col_count <- ncol(m)
      if(row_count == col_count)
      {
         total_sum <-sum(diag(m))
         total_sum
      }
      else
      {
         message ('Matrix is not square')
      }
   }
   else
   {
      message( 'Object is not a matrix')
      
   }
}
Near.Neighbor.Approach = function(x){
   x = as.data.frame(x)
   rns = x[,1]
   x = x[,-1]
   rownames(x) = rns
   x = as.matrix(x)
   
   n = nrow(x) ; k = ncol(x)
   y = x[,1]
   z = matrix(1,n,1)
   x[,1] = z
   
   prx = diag(n) - x %*% ginv(t(x)%*%x) %*% t(x)
   E = prx %*% y
   V = matrix(0,n,n)
   for (i in 1:(n-1)){
      for (j in (i+1):n){
         d = t(x[i,]-x[j,]) %*% (x[i,]-x[j,])
         if (d > 0) {
            V[i,j] = -2/d
            V[j,i] = -2/d
         }
      }
   }
   
   for (i in 1:n){
      V[i,i] = - sum(V[i,])
   }
   
   B = prx %*% V %*% prx
   
   NUM = t(E) %*% V %*% E
   DEN = t(E) %*% E/(n-k)
   STAT = NUM/DEN
   MEAN = tr(B)
   VAR = 2*((n-k)*tr(B %*% B)-(tr(B))^2)/(n-1)
   STAT = (STAT -MEAN)/sqrt(VAR)
   STAT = as.numeric(STAT)
   return(STAT)
}

t.1 = Near.Neighbor.Approach(dt1)
n = 50; k = 11
p_value = 2* pt(t.1,n-k+2)


# temporaliy elimate outliers/influential points
# temporaliy do proper transfer on X

# Variable selection

model1 = lm(M ~., data = dt1[,2:13])
k1 = ols_best_subset(model1)
plot(k1)
k1
plot(ols_stepaic_both(model1))

ols_stepwise(model1)
# choose the 9th model based on adj r^2, C(p) and AIC
dt1.a = dt1[,-c(1,6,8)]
# multicollinearity
model1.a = lm(M ~., data = dt1.a)
ols_coll_diag(model1.a)
X.corr.a = cor(dt1.a[,2:10])
corrplot(X.corr.a,method = 'circle',addCoef.col  = 'black')
# lets see if we further delete HS from dt1, how things will chagne
dt1.b = dt1[,-c(1,6,8,11)]
# multicollinearity
model1.b = lm(M ~., data = dt1.b)
ols_coll_diag(model1.b)
X.corr.b = cor(dt1.b[,2:9])
corrplot(X.corr.b,method = 'circle',addCoef.col  = 'black')
# lets see if we further delete HS from dt1, how things will chagne
dt1.c = dt1[,-c(1,3,6,7,8,10)]
# multicollinearity
model1.c = lm(M ~., data = dt1.c)
ols_coll_diag(model1.c)
X.corr.c = cor(dt1.c[,2:7])
corrplot(X.corr.c,method = 'circle',addCoef.col  = 'black')
# lets see if we further delete HS from dt1, how things will chagne
dt1.d = dt1[,-c(1,4,6,8,9)]
# multicollinearity
model1.d = lm(M ~., data = dt1.d)
ols_coll_diag(model1.d)
X.corr.d = cor(dt1.d[,2:7])
corrplot(X.corr.d,method = 'circle',addCoef.col  = 'black')

# regression and residual and influence diagonostic

summary(model1.a)
summary(model1.b)
summary(model1.c)
summary(model1.d)

ols_correlations(model.a)
ols_ovsp_plot(model.a)
http://www.rsquaredacademy.com/olsrr/articles/influence_measures.html
ols_diagnostic_panel(model.a)
ols_norm_test(model.a)
ols_corr_test(model.a)
ols_dfbetas_panel(model.a)
## Run LASSO 
# provide an arithmetic series
l1 = seq(from =-7, to =7, by = .001)
# conver to geometric series to use as user provided lambda
l2 = exp(-l1)
set.seed(42)

fit.lasso.cv = cv.glmnet(x = as.matrix(dt1[,2:12]), y =as.matrix(dt1[,1]),type.measure = 'mse', alpha = 1, family = 'gaussian', nfolds=5,  lambda = l2)
fit.ridge.cv = cv.glmnet(x = as.matrix(dt1[,2:12]), y =as.matrix(dt1[,1]),type.measure = 'mse', alpha = 0, family = 'gaussian', nfolds=5,  lambda = l2)
par(mfrow=c(1,3))
plot(fit.lasso.cv)
plot(fit.lasso.cv$glmnet.fit, "norm",   label=TRUE)
plot(fit.lasso.cv$glmnet.fit, "lambda", label=TRUE)
par(mfrow=c(1,1))
plot(fit.ridge.cv)
fit.lasso.cv$lambda.1se







#
dt2 = cbind(dt[,15],dt[,-15])
dt2 = dt2[,c(2,1,3:27)]
dt2[,2:27] = dt2[,2:27][,lapply(.SD,as.numeric)]
# standardize all X
dt2[,3:27] = dt2[,3:27][,lapply(.SD,scale)]

t.2 = Near.Neighbor.Approach(dt2)
n = 50; k = 25
p_value = 2* pt(t.2,n-k+2)

dt2.a = dt2
dt2.a[,2] = 1/dt2.a[,2]

t.2.a = Near.Neighbor.Approach(dt2.a)
n = 50; k = 25
p_value.2a = 2* pt(t.2.a,n-k+2)
colnames(dt2.a)[2] = 'inv_MA'
model2.a = lm(inv_MA ~., data = dt2.a[,2:27])

dt2.b = dt2.a
dt2.b[,6:8] = dt2.b[,6:8][,lapply(.SD,rank)]
dt2.b[,6:8] = dt2.b[,6:8]/50
dt2.c = dt2.b[,c(1:8, 12,14:16, 19:20, 23, 26:27)]
model2.c = lm(inv_MA ~., data = dt2.c[,2:ncol(dt2.c)])

time = proc.time()
k2.c = ols_best_subset(model2.c)
proc.time() - time

ols_norm_test(model2.d)


g8 = list()
for (i in 3:ncol(dt2.a)){
  temp <- ggplot(dt2.a, aes_string(x=colnames(dt2.a)[i], y = colnames(dt2.a)[2])) + labs(y ='1/MA') +
    geom_point() + theme_classic() + geom_smooth(method = 'lm')+ theme(text = element_text(size=26, face = 'bold'))
  g8[[i-2]] <- ggMarginal(temp,type = 'histogram', fill = "cyan")
}
ggarrange(g8[[1]],g8[[2]],g8[[3]],g8[[4]],g8[[5]],g8[[6]],g8[[7]],
          g8[[8]],g8[[9]],g8[[10]],g8[[11]],g8[[12]],g8[[13]],g8[[14]],
          g8[[15]],g8[[16]],g8[[17]],g8[[18]],g8[[19]],g8[[20]],
          g8[[21]],g8[[22]],g8[[23]],g8[[24]],g8[[25]],ncol = 5,nrow = 5)

model2 = lm(MA ~., data = dt2[,2:27])
k2 = ols_best_subset(model2)
plot(k2)
k2
ols_stepaic_both(model2)
ols_stepwise(model2)



Near.Neighbor.Approach(dt2)



http://www.biostat.jhsph.edu/~iruczins/teaching/jf/ch10.pdf
https://stats.stackexchange.com/questions/61217/transforming-variables-for-multiple-regression-in-r
http://nymetro.chapter.informs.org/prac_cor_pubs/01-08%20INFORMS_Jan2008.pdf
http://www.stat.columbia.edu/~martin/W2024/R10.pdf
https://stats.stackexchange.com/questions/214682/stepwise-regression-in-r-how-does-it-work
http://www.columbia.edu/~so33/SusDev/Lecture_6.pdf




## directly do v selection based on R square

