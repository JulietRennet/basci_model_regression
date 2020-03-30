##############################################
#To select a good model with only one variable check:
# 1. F-statistic and p-value: F-statistic needs to be large, depending of n. for large n, F can be slightly more that 1, else it needs to be very large.
#       The pvalue gives the significance of Fstatistic given n
# 2. R2 : the fraction of variance explained. In a simple model y ~x is the correlation between x,y. Adjusted R2 is used for polynomial functions.
# 3. RSE : RSE = sqr(RSS/(n-2)), where RSS is residual sum of squares. Is roughly how much your estimate deviates from the true regression line.
#       It is an absolute measure on how well the model fits the data.
# 4. Check the residual plot: shouldn't have any shape and the points should all be close to the 0 horizontal line. If the residuals show a 'funnel shape',
#       a transformation for the data is reccomended - log transform is the classical. Also outliers can show in the plot (leverage points)
# 5. MAE: mean absolute error - to be looked at in the test set : cross validation
# 6. RMSE: root mean squared error - to be looked at in the test set : cross validation
# 7. Bootstrapping : estimate the uncertainty of the coefficient and the model
##############################################

# In model selection you first define the model to test (linear should always be the gold standard)
# First thing always look at your data!
# Look at the F-statistic, R2 and RSE that give you a statistical measure on how well the model fit the data
# Look at the statistical test to asses if y and x have a relationship  
# Look at the diagnostic plots to identify biases and acuracy
# Look at the RMSE, R2 and MAE in the test set via cross validation
#Look at the train error and uncertainty of your coefficients with bootstrapping (at the end of the script!)


# import libraries
library(nlme) # for regression
library(ggplot2) # for plotting
library(mvtnorm)
library(boot)
library(dplyr)
library(readxl)
library(caret) #for cross validation
# Set your working directory
setwd('/home/giulia/Documents/projects/growth/model_R')


#Import data
dataset <- read_excel('phenospex.xlsx', sheet = 'DigitalBiomassvsFWvsDW') # change file name and sheet number/name
dataset <- read_excel('phenospex.xlsx') # change file name and sheet number/name

#Tranform the data
digital_biomass=as.numeric(dataset[['Digital biomass [mmÂ³]']])
FW = as.numeric(dataset[['FW (g)']])
DW = as.numeric(dataset[['DW (g)']])

df_pre <- data.frame(digital_biomass,FW,DW, stringsAsFactors = FALSE)
#Define variables
y = 'digital_biomass'
x = 'FW' #FW, DW
#prepare the dataframe
df <- df_pre[c(x,y)]
colnames(df) <- c('x', 'y')


#####################################################################
# If you want to resample your data use the following
#df = sample_frac(df, 0.2) # takes the fracion that you want to resample
####################################################################



#Investigate data
ggplot(df, aes(x = df$x, y = df$y)) + 
  geom_point() + 
  labs(x = x, y= y, title = 'scatter plot') + theme_bw()

#set polynomial formula
formula_poly = 'y~poly(x,2)'

###########################################################################################################################################
#Model fitting
# Fit model
#Linear model
linear_mod = lm(y~x, data = df) #fit y onto x

linear = fortify(linear_mod)
summary_linear = summary(linear_mod) #Gives an overview of the fit
summary_linear #Check F-statistic and p-value, tstatistic for your coefficient

#polynomial model
polynomial_mod = lm(formula_poly, data = df) #can be changed, currently set to quadratic. keep in find that going to high polynomial overfits
polynomial = fortify(polynomial_mod)
summary_polynomial = summary(polynomial_mod)
summary_polynomial #Check F-statistic and p-value, tstatistic for your coefficient


#Save residuals in data frame
df$residuals_linear <- linear$.resid
df$residuals_polynomial <- polynomial$.resid

df$st_residuals_linear <- linear$.stdresid
df$st_residuals_polynomial <- polynomial$.stdresid

df$fitted_linear <- linear$.fitted
df$fitted_ploynomial <- polynomial$.fitted

###########################################################################################################################################
#Diagnostic plotting
#Plot the data
##Linear
#Fitted line
#Plots the regression line of your model and highlights the points in the original dataframe that deviate the most from it
ggplot(df, aes(x=df$x, y=df$y)) +
  geom_point(aes(color = abs(df$residuals_linear))) + scale_color_continuous(low = 'black', high = 'blue') + #color per distance from 
  guides(color = FALSE) + #no legend for colors
  geom_segment(aes(xend = df$x, yend = df$fitted_linear), alpha = .2) +
  geom_line(aes(y=df$fitted_linear), color = 'lightgrey')+
  labs(x=x, y=y, title= ' Linear fit')  + theme_bw()

#Residuals
# Plots the residuals per fitted values and a regression line on top. Residuals should all be around mean 0 and should not have a pattern
ggplot(df, aes(x=df$fitted_linear, y=df$residuals_linear)) +
  geom_point() + 
  geom_hline(yintercept = 0, col = 'darkred', linetype = 'dashed')  +
  stat_smooth(method = 'loess') +
  labs(x='fitted values', y='residuals', title= ' Residuals vs Fitted Linear')  + theme_bw()

#QQ Plot -  to check normality assumption
#Plots the quantiles quantiles to check distribution
ggplot(df, aes(sample=st_residuals_linear))+
   stat_qq() + stat_qq_line() +
  labs(x = 'Theretical quantiles', y = 'Standardized residuals', title = 'Normal Q-Q plot linear') + theme_bw()

# Histogram of the residuals
# check if they are normally distributes
ggplot(df, aes(x = df$residuals_linear)) + 
  geom_histogram(aes(y=..density..), color='black', fill='grey') + geom_density(aes(y = ..density..), color = 'darkred') +
  labs(x = 'Residuals', y = 'density', title = 'Residuals histogram linear') + theme_bw()

#Scale location. Residuals are standardised to account for difference in variance for different observations
# spread of residuals per prediction, to check homoscedasticity (line should be somehow horizontal) and points over 2 should be investigated as possible outliers
ggplot(df, aes(x=df$fitted_ploynomial, y=sqrt(abs(df$st_residuals_linear)))) +
  geom_point()+
  stat_smooth(method='loess')+
  labs(x = 'fitted values', y = expression(sqrt("|Standardized residuals|")), title = 'Scale location plot linear') + theme_bw()

#---------------------------------------
#Polynomial
#Fitted line
#Plots the regression line of your model and highlights the points in the original dataframe that deviate the most from it
ggplot(df, aes(x=df$x, y=df$y)) + 
  geom_point(aes(color = abs(df$residuals_polynomial))) + scale_color_continuous(low = 'black', high = 'red') + #color per distance from 
  guides(color = FALSE) + #no legend for colors
  geom_segment(aes(xend = df$x, yend = df$fitted_ploynomial), alpha = .2) +
  geom_line(aes(y=df$fitted_ploynomial), color = 'lightgrey')+
  labs(x=x, y=y, title= 'Polynomial fit')  + theme_bw()

#Residuals
# Plots the residuals per fitted values and a regression line on top. Residuals should all be around mean 0 and should not have a pattern
ggplot(df, aes(x=df$fitted_ploynomial, y=df$residuals_polynomial)) +
  geom_point() + 
  geom_hline(yintercept = 0, col = 'darkred', linetype = 'dashed')  +
  stat_smooth(method = 'loess') +
  labs(x='fitted values', y='residuals', title= ' Residuals vs Fitted Polynomial')  + theme_bw()

#QQ Plot -  to check normality assumption
#Plots the quantiles quantiles to check distribution
ggplot(df, aes(sample=st_residuals_polynomial))+
  stat_qq() + stat_qq_line() +
  labs(x = 'Theretical quantiles', y = 'Standardized residuals', title = 'Normal Q-Q plot Polynomial') + theme_bw()

# Histogram of the residuals
# check if they are normally distributes
ggplot(df, aes(x = df$residuals_polynomial)) + 
  geom_histogram(aes(y=..density..), color='black', fill='grey') + geom_density(aes(y = ..density..), color = 'darkred') +
  labs(x = 'Residuals', y = 'density', title = 'Residuals histogram Polynomial') + theme_bw()

#Scale location. Residuals are standardised to account for difference in variance for different observations
# spread of residuals per prediction, to check homoscedasticity (line should be somehow horizontal) and points over 2 should be investigated as possible outliers
ggplot(df, aes(x=df$fitted_ploynomial, y=sqrt(abs(df$st_residuals_polynomial)))) +
  geom_point()+
  stat_smooth(method='loess')+
  labs(x = 'fitted values', y = expression(sqrt("|Standardized residuals|")), title = 'Scale location plot polynomial') + theme_bw()


###########################################################################################################################################
##########################################################################################################################################
#####Statistical test for model fitting between different models

#Check which model performs better. H0: models fit equally well, H1 polynomialmodel is superior
anova(linear_mod, polynomial_mod)

#Check normality distribution of residuals 
shapiro.test(df$residuals_linear)
shapiro.test(df$residuals_polynomial)

#Check values of fitting

stats = data.frame(index = c('linear', 'polynomial'),
           "Rsquared" = c(summary_linear$r.sq, summary_polynomial$r.sq),
           "adj_Rsuqared" = c( summary_linear$adj.r.squared, summary_polynomial$adj.r.squared), # included for polynomials
           "RSE"= c(summary_linear$sigma,summary_polynomial$sigma),
          "BIC" = c(BIC(linear_mod), BIC(polynomial_mod)),
          "AIC" = c(AIC(linear_mod), AIC(polynomial_mod)),

          'residuals_normality_pvalue' = c(shapiro.test(df$residuals_linear)[[2]], shapiro.test(df$residuals_polynomial)[[2]]) #high is normally distributed
           )
stats


#######################
#To check coefficient estimates and related p-value 
# It informs us on the coefficient value for the model. the pvalue from the ttest gives as the significance of the reletionship between the coefficient x and y
summary_linear$coefficients
summary_polynomial$coefficients

########################################################################################################################## 
##Cross validation and bootstrapping



# K-fold cross validation : 
#Kfold cross validation assesses the performance of the model on unseen data. It return the test error, Rsquared and MAE
# The crossvalidation iteratively divides 10 times the dataset in train (fits the model) and test set (predict with the model) and records the results
#  on which it calculates the error. This operation is repeated 3 times in the default settings, that can be changes

set.seed(42)  #Set the random seed
train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 3)  #Performs 10-fold cross validation in 3 repeats
cv_linear <- train(y~x, data = df, method = "glm", 
               trControl = train.control)

cv_poly <- train(as.formula(formula_poly), data = df, method = 'glm',
                   trControl = train.control)

cv_linear
cv_poly


# Bootstrapping :
# Randomly resample the dataset and fits the model. It provides an indication of the variance of the performance of the model
# It outputs the average model performance across the subsamples 
set.seed(42) 
train.control <- trainControl(method = "boot", number = 100) # number of bootstrapping performed
boot_linear <- train(y~x, data = df, method = "glm", 
                   trControl = train.control)

boot_poly <- train(as.formula(formula_poly), data = df, method = 'glm',
                 trControl = train.control)

#Estimate uncertainty of coefficients with bootstrapping
coef <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- glm(formula, data=d)
  return(fit$coefficients)
} 

boot_coef_linear <- boot(data = df, statistic = coef, R=1000, formula = y~x) #calculate the standard error on coeficients estimates
boot_coef_polynomial <-  boot(data = df, statistic = coef, R=1000, formula = formula_poly)

##Confidence interval for the coefficients is  coef - 2*SE(coef) -> it is the range of values that jas 95% probability of containing the true value



