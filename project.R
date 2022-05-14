install.packages("leaps")
install.packages("nclSLR", repos="http://R-Forge.R-project.org")
install.packages("mlbench")
#Load mlbench package
library(mlbench)
## Load the data
data(BreastCancer)

#Check size
dim(BreastCancer)
head(BreastCancer)
sapply(BreastCancer,class)

#choosing which columns we want to change into variables 
i = c(1:10)
# converting the factors to variables
#We need to do that because most of the statistical concepts  can only be applied 
BreastCancer[ , i] = apply(BreastCancer[ , i], 2, function(x) as.numeric(as.character(x)))
#Checking the class of each column 
sapply(BreastCancer,class)

#Now deleting the rows with NA values and creating a new data frame 
NewBreastCancer = na.omit(BreastCancer[, 2:11])
dim(NewBreastCancer)
sapply(NewBreastCancer,class)

# Explaratory analysis 
install.packages("GGally")
library(GGally)
ggpairs(Breast_Cancer_data, lower = list(combo = wrap("facethist", bindwidth = 0.5)), aes(colour = as.factor(y)))

# We can produce this advanced pairs plot which show he corelation between two variables and is
# colored based on whether the tumor is benign or malignant 

#Just by looking at the pairwise graph we see that there is a strong positive relationship 
# between cell size and cell shape so it is unlikely that we will need both predictors 
# in our model. We see that the correlation between these two predictors is 0.91 
# Also by observing the correlation matrix we observe a strong relationship between 
#cell size and cell shape with Bl chromatin and normal nucleoili, so we may not need all 
#four predictors in the model. We observe some stong relationship between some additinal
#predicotr. We will discuss these more as we go on. 

table(NewBreastCancer$Class)

#We observe that from the 683 patients 239 (35%) were found to have malignant tissue samples

apply(NewBreastCancer[,1:9], 2, mean)

# Now by looking at the predictor variables means, we see that the cell thickness is the 
# one with the highest, and mitoses is the one with the lowest. So we see that in general 
# observe more people with tissues that have higher cell size. 

# Before we start with the analysis it would be useful to standardise the predictor variables
# in order to vary on the same scale. This is important for the interpretation of our 
#results and for consistency in our analysis 

# Extract response variable
# We also need to convert the class to a numeric factor where 1 is for benign and 2 for
#malignant 
y = as.numeric(NewBreastCancer[,10]) -1
# Extract predictor variables
X1orig = NewBreastCancer[,-10]
# Standardise predictor variables
X1 = scale(X1orig)
## Combine response and standardised predictors variables in a new data frame
Breast_Cancer_data = data.frame(X1,y)
head(Breast_Cancer_data)
p = ncol(X1)
n = nrow(X1)
table(Breast_Cancer_data$y)

round(apply(Breast_Cancer_data[Breast_Cancer_data[, "y"] == 0,], 2, mean),3)
round(apply(Breast_Cancer_data[Breast_Cancer_data[, "y"] == 1,], 2, mean),3)
# We will compare the predictive performance of the different techniques using 
#10-fold cross-validation. To make the comparison fair, we will use the same set of 10 folds 
#in each case. We can sample and store the 10-folds as follows:

set.seed(1)
nfolds = 10
## Set the seed to make the analysis reproducible

## Sample fold-assignment index
fold_index = sample(nfolds, n, replace=TRUE)


# Building clasifiers 

#Of the three subset selection methods considered in lectures, best subset selection is 
#to be preferred if n > p . This is because, unlike the automated (stepwise) selection 
#methods, it considers all possible subsets of variables and so is guaranteed to find the 
#global optimum under the chosen model comparison criterion. 
#In our case we have p = 10 and n = 683 and so best subset selection can be used.

library(leaps)
library(glmnet)
library(bestglm)
## SUBSET SELECTION
glm_4 = glm(y ~ Cl.thickness + Cell.shape + Bare.nuclei +  Bl.cromatin , data = Breast_Cancer_data, family = binomial)
summary_glm_4 = summary(glm_4)
summary_glm_4$coefficients

## Apply the best subset selection method using BIC
bss_BIC_fit = bestglm(Breast_Cancer_data, family = binomial, method="exhaustive", nvmax=p) 

## Summarise the results
bss_BIC_summary = bss_BIC_fit$Subsets 
bss_BIC_summary

(best_BIC = which.min(bss_BIC_summary$BIC)-1)

# We observe that when using BIC, the 6th model is selected which includes 
#cell thickness,  marg.adhesion, bare.nuclei, Bl.cromatin,  and normal nucleoni

# Extracting the coefficients of the chosen model 
bss_BIC_fit$BestModel
bss_BIC_fit$Subsets


## Apply the best subset selection method using AIC
bss_AIC_fit = bestglm(Breast_Cancer_data, family = binomial, method="exhaustive", nvmax=p, IC = "AIC")  

bss_AIC_summary = bss_AIC_fit$Subsets
bss_AIC_summary

(best_AIC = which.min(bss_AIC_summary$AIC) -1)

# We observe that when using the AIC the 8th model is selected, which includes 
#cell thickness, cell shape, marg.adhesion, bare.nuclei, Bl.cromatin, normal nucleoni
# and mitoses


# Then compute the MSE
reg_cv = function(X1, y, fold_ind) {
  Xy = data.frame(X1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for(fold in 1:nfolds) {
    glm_fit = glm(y ~ ., data=Xy[fold_ind!=fold,], family = binomial)
    phat = predict(glm_fit, Xy[fold_ind==fold,], type = "response")
    yhat = ifelse(phat > 0.5, 1, 0) 
    yobs = y[fold_ind == fold]
    cv_errors[fold] = 1 - mean(yobs == yhat)
  }
  fold_sizes = numeric(nfolds)
  for(fold in 1:nfolds){
    fold_sizes[fold] = length(which(fold_ind==fold))
    test_error = weighted.mean(cv_errors, w=fold_sizes)
    return(test_error)
  }
}

glm_final_mse = reg_cv(X1, y, fold_index)
glm_final_mse

reg_bss_cv = function(X1, y, best_models, fold_index) {
  p = ncol(X1)
  test_errors = numeric(p)
  for(k in 1:p) {
    test_errors[k] = reg_cv(X1[,best_models[k,]], y, fold_index)
  }
  return(test_errors)
}

## Apply the function to the data
bss_mse = reg_bss_cv(X1, y, as.matrix(bss_BIC_fit$Subsets[2:10,2:10]), fold_index)
## Identify model with the lowest error
(best_cv = which.min(bss_mse))

## Create multi-panel:
par(mfrow=c(1,3))
## Produce plots, highlighting optimal value of k:
plot(1:9, bss_BIC_summary$BIC[2:10], xlab="Number of predictors", ylab="BIC", type="b")
points(best_BIC, bss_BIC_summary$BIC[best_BIC +1], col="red", pch=16)
plot(1:9, bss_AIC_summary$AIC[2:10], xlab="Number of predictors", ylab="AIC", type="b")
points(best_AIC, bss_AIC_summary$AIC[best_AIC +1], col="red", pch=16)
plot(1:9, bss_mse, xlab="Number of predictors", ylab="Test error", type="b")
points(best_cv, bss_mse[best_cv], col="red", pch=16)


glm_fit = glm(y ~ ., 
                data = Breast_Cancer_data, family = binomial)
round(glm_fit$coefficients, 3)


# Even though the different methods select different models, if we examine the plots for BIC and AIC
#we see little difference between models M5 and M8. 


## REGULARISATION METHODS 
# RIDGE REGRESSION 

## Choose grid of values for the tuning parameter
library(glmnet)
grid = 10^seq(5, -3, length=500)
## Fit a ridge regression model for each value of the tuning parameter 
ridge_fit = glmnet(X1, y, alpha=0, standardize=FALSE, lambda=grid, family = "binomial")

par(mfrow=c(1,2))
plot(ridge_fit, xvar="lambda", col=1:10, label=TRUE)

# We need to choose the appropriate tuning parameter 
# Compute 10-fold cross-validation error by usinf the same folds as in subset selection
ridge_cv_fit = cv.glmnet(X1, y, alpha=0, standardize=FALSE, lambda=grid, nfolds=nfolds, foldid=fold_index,
                         family = "binomial", type.measure = "class")
# Examine the effect of the tuning parameter on the MSE
plot(ridge_cv_fit)

#Now we can identify the optimal value for the tuning parameter
(lambda_ridge_min = ridge_cv_fit$lambda.min)

# The optimal value is lamda = 0.01535935
which_lambda_ridge = which(ridge_cv_fit$lambda == lambda_ridge_min)

## Find the parameter estimates associated with optimal value of the tuning parameter 
round(coef(ridge_fit, s=lambda_ridge_min),3)

round(summary_glm$coefficients, 3)


# We observe that most coefficients have shrunk to zero. If we compare the coefficients 
# with the ones of the full model we see that for all the variables that

# Now the corresponding cross validation error is:
(ridge_final_mse = ridge_cv_fit$cvm[which_lambda_ridge])

# Now the test error is  0.03074671
# Comparing the error with that of subset selection we see that they are almost the same 
# with this one being slightly higher. 

phat = predict(ridge_fit,s = lambda_ridge_min, newx = X1, type="response")
## Compute fitted (i.e. predicted) values:
yhat = ifelse(phat > 0.5, 1, 0)
## Calculate confusion matrix:
(confusion = table(Observed= y, Predicted=yhat))
# The training error is:
1 - mean(y == yhat)


## THE LASSO
## Choose grid of values for the tuning parameter
grid1 = 10^seq(5, -3, length=500)
## Fit a LASSO regression for each value of the tuning parameter 
lasso_fit = glmnet(X1, y, alpha=1, standardize=FALSE, lambda=grid1, family = "binomial")

## Examine the effect of the tuning parameter on the parameter estimates 
plot(lasso_fit, xvar="lambda", col=1:10, label=TRUE)

# Examining the plot we see that variable 9 is the first to
# drop out, followed by variables 5 and 9. The last variables to drop out are 2, 3 and 6 

# Compute 10-fold cross-validation error using the same folds as in ss and ridge regression
lasso_cv_fit = cv.glmnet(X1, y, alpha=1, standardize=FALSE, lambda=grid, nfolds=nfolds, foldid=fold_index,
                         family = "binomial", type.measure = "class")
plot(lasso_cv_fit)

## Identify the optimal value for the tuning parameter
(lambda_lasso_min = lasso_cv_fit$lambda.min) 
which_lambda_lasso = which(lasso_cv_fit$lambda == lambda_lasso_min)
## Find the parameter estimates associated with optimal value of the tuning parameter
round(coef(lasso_fit, s=lambda_lasso_min), 3)

# We see that the coefficients for id has shrunk to zero. However for all
# the other coefficients we do not observe any significant reduction 

# Now the corresponding cross validation error is:
(lasso_final_mse = lasso_cv_fit$cvm[which_lambda_lasso])

# The test error now is 0.03367496 



## DISCRIMINANT ANALYSIS 
# LDA 
library("plyr")
library(nclSLR)
library(tidyr)
library(combinat)

# Create a vector that contain the names of all the predictor variables 

my_vec = c("Cl.thickness",  "Cell.size", "Cell.shape", "Marg.adhesion",
           "Epith.c.size", "Bare.nuclei", "Bl.cromatin", "Normal.nucleoli","Mitoses")
# Now compute all the possible subsets of the variables 
my_combi1 = unlist(lapply(1:length(my_vec),  combinat::combn,  x = my_vec, simplify = FALSE), recursive = FALSE)
my_combi1


install.packages("MASS")
library(MASS)

# Create the function that performs cross validation on lda 
reg_cv_lda = function(x1, y, fold_ind){
  Xy = data.frame(x1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for (fold in 1:nfolds) {
    tmp_fit = lda(y~., data = Xy[fold_ind!=fold,])
    phat = predict(tmp_fit, Xy[fold_ind == fold,])
    yhat = phat$class
    yobs = y[fold_ind==fold]
    cv_errors[fold] = 1 - mean(yobs == yhat)
  }
  fold_sizes = numeric(nfolds)
  for (fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error_lda = weighted.mean(cv_errors, w=fold_sizes)
}

# Now create a for loop that runs the function for every possible subset and finds the cv error 
errors_lda = numeric(length(my_combi1))
for (k in 1:length(my_combi1)){
  df = as.data.frame(Breast_Cancer_data[, unlist(my_combi1[k])])
  errors_lda[k] = reg_cv_lda(df,y,fold_index)
}
errors_lda

# Find which subset had the lowest cv error 
which.min(errors_lda)
# This is the 413 subset
# Find the specific cv errors 
errors_lda[413]
# The test error for this model is 0.03660322
# Now find the predictor variables that this model contains
my_combi1[413]
# So this model contains "Cl.thickness" "Cell.size"    "Epith.c.size" "Bare.nuclei"  "Bl.cromatin"  "Mitoses" 

# Now we can perform LDA and find the group means 
llda = lda(y~ Cl.thickness+ Cell.size + Epith.c.size + Bare.nuclei + Bl.cromatin + Mitoses, 
           data = Breast_Cancer_data)
llda$means
llda

# QDA 
# Now we can perform qda using the same way as before
# Create the function that performs cross validation on qda 
reg_cv_qda = function(x1, y, fold_ind){
  Xy = data.frame(x1, y=y)
  nfolds = max(fold_ind)
  if(!all.equal(sort(unique(fold_ind)), 1:nfolds)) stop("Invalid fold partition.")
  cv_errors = numeric(nfolds)
  for (fold in 1:nfolds) {
    tmp_fit = qda(y~., data = Xy[fold_ind!=fold,])
    phat = predict(tmp_fit, Xy[fold_ind == fold,])
    yhat = phat$class
    yobs = y[fold_ind==fold]
    cv_errors[fold] = 1 - mean(yobs == yhat)
  }
  fold_sizes = numeric(nfolds)
  for (fold in 1:nfolds) fold_sizes[fold] = length(which(fold_ind==fold))
  test_error_lda = weighted.mean(cv_errors, w=fold_sizes)
}

# Now create a for loop that runs the function for every possible subset and finds the cv error 
errors_lda = numeric(length(my_combi1))
for (k in 1:length(my_combi1)){
  df = as.data.frame(Breast_Cancer_data[, unlist(my_combi1[k])])
  errors_qda[k] = reg_cv_qda(df,y,fold_index)
}
errors_qda

# Find which subset had the lowest cv error 
which.min(errors_qda)
# This is the 166 subset
# Find the specific cv errors 
errors_qda[166]
# The test error for this model is 0.03513909
# Now find the predictor variables that this model contains
my_combi1[166]
# This model contains "Cl.thickness"  "Marg.adhesion" "Epith.c.size"  "Bare.nuclei"  

# Now we can perform qda on this model 
quadraticDA = qda(y~ Cl.thickness+ Marg.adhesion + Epith.c.size  + Bare.nuclei, 
               data = Breast_Cancer_data)
quadraticDA$means

