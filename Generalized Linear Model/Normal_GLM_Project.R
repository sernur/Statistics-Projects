load("Normal_GLM.rdata")
source("Normal_GLM.r")
library(carData)
library(car)

# (a)
xsquared = x^2
fit = lm(y~x+xsquared)
summary(fit)
# plot(fit)

# Test whether the quadratic term is significant... p-value is in the printout though
linearHypothesis(fit, "xsquared = 0")

# fitH0 = lm(y~x)
# SSres = sum(fit$residuals^2)
# SSresH0 = sum(fitH0$residuals^2)
# 
# F = ((SSresH0 - SSres) / (SSres)) * ((50-3)/(1))
# 1-pf(F, 1, 50-3)

# Plot the residuals and decide if the data is homoscedastic.... Possibly heteroscedastic...
plot(fit$residuals)



# (b) calculating the log-likelihood 

beta1 = fit$coef / var(y)  # eta_1 is mu_1 / sigma_i^2 and in the earlier model the sigma is constant (sample variance y)
beta2 = c(1 / var(y), 0, 0)  # eta_2 is 1/sigma_i^2 and in the earlier model sigma is constant and doesn't depend on x at all
ones = rep(1, length(x))
X1 = cbind(ones, x, x^2) 
X2 = cbind(ones, x, x^2)
lik(beta1, beta2, y, X1, X2)

# gradient 
Dlik(beta1, beta2, y, X1, X2)  # note the gradient for eta_1 is very small, for eta_2 the gradient is way off. The way to improve is to change the variance estimates to depend on x.  We are currently at a maximum where the sigma doesn't depend on x


# two explanations:  The model in part (a) is a submodel where beta_2,1 = beta_2,2 = 0 and have the tranformations above for the betas.
# Gradient explanation: the gradient for eta_1 is very small, for eta_2 the gradient is way off. The way to improve is to change the variance estimates to depend on x.  We are currently at a maximum where the sigma doesn't depend on x


# (c) Optimizing

# Optimizing with Newton-Raphson
NewtonRaphsonOptim = function(currentBeta) {
  while (TRUE) {
    l = lik(currentBeta[1:3], currentBeta[4:6], y, X1, X2)
    Dl = Dlik(currentBeta[1:3], currentBeta[4:6], y, X1, X2)
    D2l = D2lik(currentBeta[1:3], currentBeta[4:6], y, X1, X2)
    
    if (sum(Dl^2) < 10^(-8)) {
      return(currentBeta)
    }
    
    hm = -solve(D2l) %*% Dl
    newBeta = currentBeta + hm
    newl = lik(newBeta[1:3], newBeta[4:6], y, X1, X2)
    
    while (is.nan(newl) || newl <= l) {
      hm = hm / 2
      newBeta = currentBeta + hm
      newl = lik(newBeta[1:3], newBeta[4:6], y, X1, X2)
    }
    
    currentBeta = newBeta
  }
}

betaMLE = NewtonRaphsonOptim(c(beta1, beta2))
betaMLE

# Optimizing with optim()
optimized = optim(c(beta1, beta2), 
      function(betaC) {lik(betaC[1:3], betaC[4:6], y, X1, X2)}, 
      gr=function(betaC) {Dlik(betaC[1:3], betaC[4:6], y, X1, X2)},
      method="BFGS", # This is a "quasi-Newton method (variable metric algorithm)" which takes advantage of the gradient 
      control=list(fnscale=-1))  # Maximize rather than minimize
# NOTE: I think R essentially filtered out the inappropriate beta2 values (due to NaN from the log(..)), which is why it throws these warnings
# optimized$par[4] + x*optimized$par[5] + x^2*optimized$par[6]  # The estimated beta2 makes sense though
optimized

# Plot of data vs predictions from standard linear model
plot(x, y, col= "black", ylim=c(0,16))
predictedFromSLM = fit$coef[1] + fit$coef[2]*x + fit$coef[3]*(x^2) 
points(x, predictedFromSLM, col="#CC79A7") # pink

estimatedSigmaSq = 1/(betaMLE[4] + betaMLE[5]*x + betaMLE[6]*(x^2))
predictedFromGLM = (betaMLE[1] + betaMLE[2]*x + betaMLE[3]*(x^2)) * estimatedSigmaSq
points(x, predictedFromGLM, col="#0072B2") # blue


residualsGLM = predictedFromGLM - y
plot(x, residualsGLM, ylim=c(-1.5,1.5))
#lines(x, residualsGLM+sd(residualsGLM))
#lines(x, residualsGLM-sd(residualsGLM))
#lines(x, residualsGLM+sqrt(estimatedSigmaSq), col="#0072B2")
#lines(x, residualsGLM-sqrt(estimatedSigmaSq), col="#0072B2")
lines(x, sqrt(estimatedSigmaSq), col="#0072B2")
lines(x, -sqrt(estimatedSigmaSq), col="#0072B2")

sqrt(estimatedSigmaSq)+residualsGLM

residualsSLM = predictedFromSLM - y
plot(x, residualsSLM, ylim=c(-5,4))
#lines(x, residualsSLM+1.08)
#lines(x, residualsSLM-1.08)
lines(x, rep(1.08, length(x)), col="#0072B2")
lines(x, rep(-1.08, length(x)), col="#0072B2")


# (d) 
# test for quadratic term in eta_1 
NROptimNoQuadEta1 = function(currentBeta) {
  finished = FALSE
  while (!finished) {
    l = lik(currentBeta[1:2], currentBeta[3:5], y, X1[,1:2], X2)
    Dl = Dlik(currentBeta[1:2], currentBeta[3:5], y, X1[,1:2], X2)
    D2l = D2lik(currentBeta[1:2], currentBeta[3:5], y, X1[,1:2], X2)
    
    if (sum(Dl^2) < 10^(-8)) {
      return(currentBeta)
    }
    
    hm = -solve(D2l) %*% Dl
    newBeta = currentBeta + hm
    newl = lik(newBeta[1:2], newBeta[3:5], y, X1[,1:2], X2)
    
    while (is.nan(newl) || newl <= l) {
      hm = hm / 2
      newBeta = currentBeta + hm
      newl = lik(newBeta[1:2], newBeta[3:5], y, X1[,1:2], X2)
    }
    
    currentBeta = newBeta
  }
  return(currentBeta)
}

H0betaMLE = NROptimNoQuadEta1(c(betaMLE[1:2], betaMLE[4:6]))
H0betaMLE

dev1 = -2*lik(betaMLE[1:3], betaMLE[4:6], y, X1, X2)
dev0 = -2*lik(H0betaMLE[1:2], H0betaMLE[3:5], y, X1[,1:2], X2)

1-pchisq(dev0 - dev1, df=1)
dev0-dev1

# test for quadratic term in eta_2 
NROptimNoQuadEta2 = function(currentBeta) {
  finished = FALSE
  while (!finished) {
    l = lik(currentBeta[1:3], currentBeta[4:5], y, X1, X2[,1:2])
    Dl = Dlik(currentBeta[1:3], currentBeta[4:5], y, X1, X2[,1:2])
    D2l = D2lik(currentBeta[1:3], currentBeta[4:5], y, X1, X2[,1:2])
    
    if (sum(Dl^2) < 10^(-8)) {
      return(currentBeta)
    }
    
    hm = -solve(D2l) %*% Dl
    newBeta = currentBeta + hm
    newl = lik(newBeta[1:3], newBeta[4:5], y, X1, X2[,1:2])
    
    while (is.nan(newl) || newl <= l) {
      hm = hm / 2
      newBeta = currentBeta + hm
      newl = lik(newBeta[1:3], newBeta[4:5], y, X1, X2[,1:2])
    }
    
    currentBeta = newBeta
  }
  return(currentBeta)
}

H0betaMLE = NROptimNoQuadEta2(c(betaMLE[1:3], betaMLE[4:5]))
H0betaMLE

dev1 = -2*lik(betaMLE[1:3], betaMLE[4:6], y, X1, X2)
dev0 = -2*lik(H0betaMLE[1:3], H0betaMLE[4:5], y, X1, X2[,1:2])

1-pchisq(dev0 - dev1, df=1)



# Test for GLM model vs SLM
dev1 = -2*lik(betaMLE[1:3], betaMLE[4:6], y, X1, X2)
dev0 = -2*lik(beta1, beta2, y, X1, X2)

1-pchisq(dev0 - dev1, df=2)



# (e)
predictionSimulation = function(xVal, numTrials) {
  mseSLMmean = 0
  mseGLMmean = 0
  
  mseSLMsd = 0
  mseGLMsd = 0
  
  # The true values of the mean and variance from the model in (c)
  trueVariances = estimatedSigmaSq
  trueMeans = predictedFromGLM
  # The true values for the xValue from the model in (c)
  modelVar = 1/(betaMLE[4] + betaMLE[5]*xVal + betaMLE[6]*(xVal^2))
  modelSd = sqrt(modelVar)
  modelMean = (betaMLE[1] + betaMLE[2]*xVal + betaMLE[3]*(xVal^2)) * modelVar
  
  for (trial in 1:numTrials) {
    simulatedY = rnorm(50, mean=trueMeans, sd=sqrt(trueVariances))
    newSLMfit = lm(simulatedY~x+xsquared)
    newSLMPrediction = newSLMfit$coef[[1]] + newSLMfit$coef[[2]]*xVal + newSLMfit$coef[[3]]*xVal^2
    newSLMsd = sd(newSLMfit$residuals)
    
    newBeta1 = newSLMfit$coef / var(simulatedY)  
    newBeta2 = c(1 / var(simulatedY), 0, 0)  
    newBetaMLE = NewtonRaphsonOptim(c(newBeta1, newBeta2))
    
    newGLMEstVariance = 1/(newBetaMLE[4] + newBetaMLE[5]*xVal + newBetaMLE[6]*(xVal^2))
    newPredictedY = (newBetaMLE[1] + newBetaMLE[2]*xVal 
                     + newBetaMLE[3]*(xVal^2)) * newGLMEstVariance
    newGLMsd = sqrt(newGLMEstVariance)
    
    mseSLMmean = mseSLMmean + (newSLMPrediction - modelMean)^2
    mseGLMmean = mseGLMmean + (newPredictedY - modelMean)^2
    
    mseSLMsd = mseSLMsd + (newSLMsd - modelSd)^2
    mseGLMsd = mseGLMsd + (newGLMsd - modelSd)^2
  }
  
  mseSLMmean = (1/numTrials) * mseSLMmean
  mseGLMmean = (1/numTrials) * mseGLMmean
  
  mseSLMsd = (1/numTrials) * mseSLMsd
  mseGLMsd = (1/numTrials) * mseGLMsd
  
  return(c(mseSLMmean, mseGLMmean, mseSLMsd, mseGLMsd))
}
predictionSimulation(0.75, 100)
