library(ggplot2)
library(lmtest)
library(MASS)
library(ResourceSelection)
#install.packages('gam')
library(gam)
setwd('/Users/melisa/CS7200_StatMethods/final_project/cs7200-stats-final-project')
mutations = read.csv("Data/cleaned_data.csv", header=TRUE)

# -------------------Check for Multicollinearity---------------------
# Create pairwise contingency tables for predictors to check for multicollinearity
hydro_vol = table(hydrophobicity=mutations$hydrophobicity, volume=mutations$volume) 
hydro_vol_chi = chisq.test(hydro_vol) 
hydro_hDonor = table(hydrophobicity=mutations$hydrophobicity, hDonor=mutations$hDonorAcceptor)
hydro_hDonor_chi = chisq.test(hydro_hDonor)
vol_hDonor = table(volume=mutations$volume, hDonor=mutations$hDonorAcceptor)
vol_hDonor_chi = chisq.test(vol_hDonor)
hydro_vol
hydro_vol_chi
hydro_hDonor
hydro_hDonor_chi
vol_hDonor
vol_hDonor_chi

# Run simple additive logistic regression model with volume and hydrogen donor status
# Yi = B + B1V + B2H
mutations_small <- data.frame(vol=mutations$volume, hDonor=mutations$hDonorAcceptor, disease=mutations$DiseaseOutcome)
fit_V_H <- glm(DiseaseOutcome ~ volume + hDonorAcceptor, family=binomial, data=mutations) 
summary(fit_V_H)
# confidence intervals on parameter estimates
confint(fit_V_H)

# -------------------Fitted Model from Contingency Table---------------------
# Functions to compare log-likelihood of the fitted model and the actual model 
getSampleProportion <- function(no, yes) { return(round(yes/ (no + yes), 3)) }
getPointEstimate <- function(beta, se) { return(exp(beta) / (1+exp(beta))) }
getSeLower <- function(beta, se) { return(exp(beta-se) / (1+exp(beta-se))) }
getSeUpper <- function(beta, se) { return(exp(beta+se) / (1+exp(beta+se))) }
getOdds <- function(beta) {return (exp(beta)) }
getL <- function(beta, se) {return (exp(beta - 1.96*se)) }
getU <- function(beta, se) { return (exp(beta + 1.96*se)) }

getEstimatedEffects <- function(df) {
  for (i in 1:nrow(df)) {
    df$SProp[i] <- getSampleProportion(df$No[i], df$Yes[i])
    df$PointEst[i] <- getPointEstimate(df$Beta[i], df$SE[i])
    df$SeLower[i] <- getSeLower(df$Beta[i], df$SE[i])
    df$SeUpper[i] <- getSeUpper(df$Beta[i], df$SE[i])
    df$Odds[i] <- getOdds(df$Beta[i])
    df$WaldCIEffectL[i] <- getL(df$Beta[i], df$SE[i])
    df$WaldCIEffectU[i] <- getU(df$Beta[i], df$SE[i])
  }
  return(df)
}

# Determine fitted values
ctable <- ftable(mutations_small)
ctable
summary(fit_V_H)
b0 <- fit_V_H$coefficients[1]
b1 <- fit_V_H$coefficients[2]
b2 <- fit_V_H$coefficients[3]
se0 <- 0.1700
se1 <- 0.1664
se2 <- 0.1406
V_H_Df <- data.frame(Volume=c(0, 0, 1, 1), 
                     HDonor=c(0, 1, 0, 1),
                     No=ctable[,1],
                     Yes=ctable[,2],
                     Beta=c(b0, b0+b2, b0+b1, b0+b1+b2),
                     SE=c(se0, max(se0, se2), max(se0, se1), max(se0, se1, se2))
                     )
V_H_Df <- getEstimatedEffects(V_H_Df)
V_H_Df
# plot fitted values
ggplot(data=V_H_Df, aes(x=factor(Volume), y=value, group=factor(HDonor), color=factor(HDonor)) ) + 
  geom_line(aes(y=PointEst)) +
  ylim(0.5, 1.0) +
  geom_errorbar(aes(y=PointEst, ymin=SeLower, ymax=SeUpper), width=.05) +
  geom_point(aes(y=SProp)) +
  labs(color="Hydrogen Donor Change", title="Model without Interaction", x="Volume Change", y="Probability of Disease (+-SE)")



# -------------------Assess Model Fit of Fitted Model---------------------
V_H_Df
total <- V_H_Df$Yes + V_H_Df$No
total
estimatedYes <- total*V_H_Df$PointEst
estimatedNo <- total*(1 - V_H_Df$PointEst)
estimated <- c(estimatedYes, estimatedNo)
observed<- c(V_H_Df$Yes, V_H_Df$No)
estimated
observed

pearson_stat <- sum( (observed-estimated)^2 / estimated )
# df = number of settings of x (i.e. params is saturated model) - number of params in model
pearson_p <- pchisq(pearson_stat, 4-3)
c(pearson_stat, pearson_p)

# Run G-Squared likelihood ratio aginst fitted values
# odds of disease for no change in volume and no change in H Donor Status
G_squared <- 2 * sum(observed * log(observed/estimated))
G_squared_p <- pchisq(G_squared, 4-3)
c(G_squared, G_squared_p)

# H-L Goodness of Fit (not relevant)
hoslem.test(mutations_small$disease, fitted(fit_V_H))


# Wald Test for each B_k, H0: B_k = 0
# Has Chi-Squared distribution with df=1
getWaldSummary <- function (b, se) {
  waldStat <- (b/se)^2
  pWald <- 1 - pchisq(waldStat, 1)
  return( c(waldStat, pWald) )
}
getWaldSummary(b0, se0)
getWaldSummary(b1, se1)
getWaldSummary(b2, se2)



# -------------------Goodness of Fit Compared to Other Models---------------------

# Run logistic regression models for single parmeters and 2 parameters with interaction
# Yi = B0 + B1V
fit_V <- glm(disease ~ vol, family=binomial, data=mutations_small)
# Yi = B + B1H
fit_H <- glm(disease ~ hDonor, family=binomial, data=mutations_small)
# Yi = B + B1V + B2H + B3VH
fit_V_H_VH <- glm(disease ~ vol + hDonor + vol*hDonor, family=binomial, data=mutations_small) 
summary(fit_V)
# confidence intervals on parameter estimates
confint(fit_V)
summary(fit_H)
# confidence intervals on parameter estimates
confint(fit_H)
summary(fit_V_H)
# confidence intervals on parameter estimates
confint(fit_V_H)
summary(fit_V_H_VH) #lowest AIC
# confidence intervals on parameter estimates
confint(fit_V_H_VH)


# Likelihood Ratio Tests (against the additive multiple logistic regression model)
lr1 <- lrtest(fit_V_H, fit_H)
lr2 <- lrtest(fit_V_H, fit_V) 
lr3 <- lrtest(fit_V_H, fit_V_H_VH) 
lr1
lr2
lr3



# Hosmer-Lemeshow
# H-L test of lack of fit
hoslem.test(mutations_small$disease, fitted(fit_V))
hoslem.test(mutations_small$disease, fitted(fit_H))
hoslem.test(mutations_small$disease, fitted(fit_V_H))
hoslem.test(mutations_small$disease, fitted(fit_V_H_VH))

# Compare residual deviances to chisq(1)
# compare the interaction model and no interaction model
pchisq(deviance(fit_V_H_VH)-deviance(fit_V_H), 1, lower=F)
pchisq(deviance(fit_V_H_VH)-deviance(fit_V_H), 1, lower=F)


# Pearson residuals vs predicted response
plot( residuals(fit_V_H, type="pearson") ~ predict(fit_V_H_VH, type="response"), 
      xlab=expression(hat(pi)), ylab="Pearson Residual")

# Pearson residuals vs predicted link    
plot(residuals(fit_V_H, type="pearson") ~ predict(fit_V_H_VH,type="link"), 
     xlab=expression(hat(eta)), ylab="Pearson Residual")

# Deviance residuals vs predicted response
plot( residuals(fit_V_H, type="deviance") ~ predict(fit_V_H_VH, type="response"), 
      xlab=expression(hat(pi)), ylab="Deviance Residual")

# Studentized residuals vs predicted response
plot( rstudent(fit_V_H) ~ predict(fit_V_H_VH, type="response"), 
      xlab=expression(hat(pi)), ylab="Studentized Residual")

# Cooks distance
plot(cooks.distance(fit_V_H) ~ predict(fit_V_H_VH,type="response"), 
     xlab=expression(hat(pi)), ylab="Cooks distance")



# ---------------------Evidence of mild overdispersion------------------
# Estimate the overdispersion parameter
est.phi <- function(glmobj) { 
  sum( residuals(glmobj, type="pearson")^2 / df.residual(glmobj) )
}
est.phi(fit_V_H) # if >1 evidence of overdispersion
est.phi(fit_V_H_VH) # if >1 evidence of overdispersion


# refit the interaction model while allowing for overdispersion
fit_overdispersion <- glm(disease ~ vol + hDonor + vol*hDonor, family=quasibinomial, data=mutations_small)
summary(fit_overdispersion)
summary(fit_V_H_VH)
# equivalent (more computationally efficient): update the previous model
summary( update(fit_V_H_VH, family="quasibinomial") )

# examine the resulting change in inference
# newOrings.predictLink1 <- predict(fit1, newdata=newOrings, se.fit=T, type="link")
mutations.predictLink1 <- predict(fit_overdispersion, newdata=mutations_small, se.fit=T, type="link")
L1 <- mutations.predictLink1$fit - qnorm(1-0.05/(2*10))*mutations.predictLink1$se
U1 <- mutations.predictLink1$fit + qnorm(1-0.05/(2*10))*mutations.predictLink1$se
plot(disease ~ vol + hDonor, mutations_small, xlab="volume and hdonor",ylab="Disease", pch=16)
lines(mutations$vol, 1/(1+exp(-L1)), lty=2, col="red")
lines(mutations$vol, 1/(1+exp(-U1)), lty=2, col="red")
legend("topright", lty=2, col=c("blue","red"), c("no overdispersion", "with overdispersion"))


### TAKEN FROM PROF VITEK'S EXAMPLE CODE
# --------------Predictive ability on the training set-------------------------
library(ROCR)
# select training set (as example, here only use 1/4 of the data to build the model)
train <- sample(x=1:nrow(mutations_small), size=nrow(mutations_small)/3)

# fit the full model on the training dataset
fit.train <- glm(disease ~vol + hDonor + vol*hDonor, family=binomial, data=mutations_small[train,])




summary(fit.train)
summary(fit_V_H_VH)

# calculate predicted probabilities on the same training set
scores <- predict(fit.train, newdata=mutations_small[train,], type="response")

# compare predicted probabilities to labels, for varying probability cutoffs
pred <- prediction(scores, labels=mutations_small[train,]$disease )
perf <- performance(pred, "tpr", "fpr")

# plot the ROC curve
plot(perf, colorize=F, main="In-sample ROC curve")

# print out the area under the curve
unlist(attributes(performance(pred, "auc"))$y.values)


# --------------Evalate the predictive ability on the validation set------------
# make prediction on the validation dataset
scores <- predict(fit.train, newdata=mutations_small[-train,], type="response")
pred <- prediction( scores, labels=mutations_small[-train,]$disease )
perf <- performance(pred, "tpr", "fpr")

# overlay the line for the ROC curve
plot(perf, colorize=T, add=TRUE)

# print out the area under the curve
unlist(attributes(performance(pred, "auc"))$y.values)








# --------------Run Model with Pocket Parameters------------
# Use pocket columns: 'Pocket_A', 'Pocket_G', 'Pocket_GearyAuto_AvFlexibility30'
# The data with pockets has a smaller n value

pockets = read.csv("Data/Sequential_Merged_With_Pockets.csv", header=TRUE)[ ,c('FromAA', 'ToAA', 'Pocket_A', 'Pocket_G', 'Pocket_GearyAuto_AvFlexibility30', 'Class',	'X2',	'X3',	'X2X3')]
# filter out any unknown amino acids
pockets <- pockets[(pockets$FromAA != 'X') & (pockets$ToAA != 'X'), ]
dim(pockets) #545 9
head(pockets)

# rename columns
#names(pockets) <- c()

# group the data by categories
range(pockets$Pocket_A)
range(pockets$Pocket_G)
range(pockets$Pocket_GearyAuto_AvFlexibility30)
# group every 2
pockets$pocketA = pockets$Pocket_A%/%2 * 2
pockets$pocketG = pockets$Pocket_G%/%2 * 2
# group every .2
pockets$pocketGeary = ((round(pockets$Pocket_GearyAuto_AvFlexibility30, 1) * 10) %/% 2 * 2) / 10

# look at number of bins and count per bin
table(data.frame(pockets$pocketGeary, pockets$Class))
tableA <- table(data.frame(pockets$pocketA, pockets$Class))
table(data.frame(pockets$pocketG, pockets$Class))

# plot log odds against parameter values should yield a straight line

# X’β = β0 + β1Xi1 + β2Xi2 + β3Xi1 Xi2 + β4Xi3 + β5Xi4 + β6Xi5
fit.pocket <- glm(Class ~ X2 + X3 + X2*X3 + pocketA + pocketG + pocketGeary, family=binomial, data=pockets) 
summary(fit.pocket)

# Wald CI for the parameters, 
confint(fit.pocket)

# using residuals plot, 


# Pearson residuals vs predicted response
plot( residuals(fit.pocket, type="pearson") ~ predict(fit.pocket, type="response"), 
      xlab=expression(hat(pi)), ylab="Pearson Residual")

# Pearson residuals vs predicted link    
plot(residuals(fit.pocket, type="pearson") ~ predict(fit.pocket,type="link"), 
     xlab=expression(hat(eta)), ylab="Pearson Residual")

# Deviance residuals vs predicted response
plot( residuals(fit.pocket, type="deviance") ~ predict(fit.pocket, type="response"), 
      xlab=expression(hat(pi)), ylab="Deviance Residual")





# 
# # Pearson Chi-squared test for goodness of fit 
# 
# # log likelihood ratio tests
# summary(fit_V_H_VH)
# fit.simple <- glm(Class ~ X2 + X3 + X2*X3, family=binomial, data=pockets) 
# summary(fit.simple)
# lrtest(fit.simple, fit.pocket) 
# 
# 


