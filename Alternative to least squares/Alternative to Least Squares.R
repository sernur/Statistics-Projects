### load the data
setwd("/Users/sergazy/Downloads/Fall 2018 courses and all files/Eric's class/08:30:2018 Alternative estimators")
load("dataL1L2.Rdata")
ls()

#a
fit=lm(Y~x1+x2) # in the question Y is dependent  variable is given
summary(fit)
resid<-fit$res
hist(resid)
qqnorm(resid)
u<-qqnorm(resid, plot= FALSE)
qqline(resid) # it looks a bit of.

#b
L1 = function(b,x1,x2,Y){
  sum(abs(Y-b[1]-b[2]*x1-b[3]*x2))
}
par0=fit$coef #Estimator comes from least square
par=optim(par0,L1,gr=NULL,x1,x2,Y) # estimator comes from L1 
par
#now compare with par and par0

# I am doing something here 
X=cbind(1,x1,x2)
dim(X)
X %*% as.matrix(par$par)
resid2 <- Y - X %*% as.matrix(par$par)
resid2 - resid
hist(resid2-resid)
mean(resid2 - resid)
#this is over I don't think the above is necessary

#c
X=cbind(1,x1,x2)
M=1000
b1sim=matrix(0,M,3)
b2sim=matrix(0,M,3)
error1=0
error2=0
S=sqrt(sum(fit$res^2)/47)
for (i in 1:M){
  Usim=rnorm(50)*S
  Ysim= X%*%fit$coef + Usim
  fitsim=lm(Ysim~x1+x2)
  b2sim[i,]=fitsim$coef
  parsim=optim(fitsim$coef,L1,gr=NULL,x1,x2,Ysim)
  b1sim[i,]=parsim$par
  error1=error1+sum((b1sim[i,]-fit$coef)^2)
  error2=error2+sum((b2sim[i,]-fit$coef)^2)
}
average_error1=error1/1000
average_error2=error2/1000

#d
#bla bla # Prove it 

#e
X=cbind(1,x1,x2)
M=1000
b3sim=matrix(0,M,3)
b4sim=matrix(0,M,3)
error3=0
error4=0
S=sqrt(sum(fit$res^2)/47)
for (i in 1:M){
  Usim=rt(50,3)*S
  Ysim= X%*%fit$coef + Usim
  fitsim=lm(Ysim~x1+x2)
  b4sim[i,]=fitsim$coef
  parsim=optim(fitsim$coef,L1,gr=NULL,x1,x2,Ysim)
  b3sim[i,]=parsim$par
  error4=error4+sum((b4sim[i,]-fit$coef)^2)
  error3=error3+sum((b3sim[i,]-fit$coef)^2)
}
average_error4=error4/1000 #this is L2
average_error3=error3/1000 #this is L1, this looks better which is very strange.

#f
X=cbind(1,x1,x2)
M=1000
b5sim=matrix(0,M,3)
b6sim=matrix(0,M,3)
count=0
beta_f=c(2,-1,1)
S=sqrt(sum(fit$res^2)/47)
for (i in 1:M){
  Usim=rt(50,3)*S
  Ysim= X%*%beta_f + Usim
  fitsim=lm(Ysim~x1+x2)
  b6sim[i,]=fitsim$coef #this is L2
  parsim=optim(fitsim$coef,L1,gr=NULL,x1,x2,Ysim)
  b5sim[i,]=parsim$par #this is L1
  
  S_f=vcov(fitsim)
  coef_f=fitsim$coef
  con=c(0,0.7,0)
  t_f=(coef_f[2]+1)/sqrt(t(con)%*%S_f%*%con)
  p_f=2*(1-pt(t_f,50))
  if (p_f < 0.05) {
    count=count+1
  }
  else{
    count=count+0
  }
  
}
print(count)




