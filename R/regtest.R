regtest <-
function(data, model=1, start=c(a=1,b=1,c=1,d=1,e=1)){
names(data)=c("x","treatment","y") 
data$x=as.numeric(data$x)  
data$treatment=as.factor(data$treatment) 
data$y=as.numeric(data$y)  
s=split(data,data[,2])
l <- combn(s,2)
f11o=function(x){t=rbind(l[,x][[1]],l[,x][[2]]); return(t)}
ty=nlevels(data$treatment)
i=ifelse(ty==2,1,nlevels(data$treatment))
i=1:i
r1=lapply(i,f11o)
f22=function(iiii){rr=data.frame(r1[[iiii]][,1],r1[[1]][,2],r1[[iiii]][,3]);names(rr)=c("x","treatment","y");return(rr)}
iiii=1:length(r1)
r2=lapply(iiii, f22)
    
# linear
f1=function(data){
# linear
ml1=nls(y~a+b*x, data=data, start=c(a=1,b=1))
# a
ml2=nls(y~a[treatment]+b*x, data=data, start=list(a=c(1,1),b=1))
# b
ml3=nls(y~a+b[treatment]*x, data=data, start=list(a=c(1),b=c(1,1)))
# model
ml4=nls(y~a[treatment]+b[treatment]*x, data=data, start=list(a=c(1,1),b=c(1,1)))
s1l=summary(ml1)[[10]]
s2l=summary(ml2)[[10]]
s3l=summary(ml3)[[10]]
s4l=summary(ml4)[[10]]
a1l=anova(ml1,ml2)
a2l=anova(ml1,ml3)
a3l=anova(ml1,ml4)
AIC=AIC(ml1,ml2,ml3,ml4); AIC=AIC[,-1]
BIC=BIC(ml1,ml2,ml3,ml4); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2")
l1=list(s1l,s2l,s3l,s4l,a1l,a2l,a3l, par)
names(l1)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2","test a parameter", "test b parameter","test model", "AIC and BIC") 
return(l1)
    }

f2=function(data){     
# quadratic
mq1=nls(y~a+b*x+c*x^2, data=data, start=c(a=1,b=1,c=1))
# a
mq2=nls(y~a[treatment]+b*x+c*x^2, data=data, start=list(a=c(1,1),b=c(1), c=c(1)))
# b
mq3=nls(y~a+b[treatment]*x+c*x^2, data=data, start=list(a=c(1),b=c(1,1), c=c(1)))
# c
mq4=nls(y~a+b*x+c[treatment]*x^2, data=data, start=list(a=c(1),b=c(1), c=c(1,1)))
# model
mq5=nls(y~a[treatment]+b[treatment]*x+c[treatment]*x^2, data=data, start=list(a=c(1,1),b=c(1,1), c=c(1,1)))
s1q=summary(mq1)[[10]]
s2q=summary(mq2)[[10]]
s3q=summary(mq3)[[10]]
s4q=summary(mq4)[[10]]
s5q=summary(mq5)[[10]]
a1q=anova(mq1,mq2)
a2q=anova(mq1,mq3)
a3q=anova(mq1,mq4)
a4q=anova(mq1,mq5)
AIC=AIC(mq1,mq2,mq3,mq4, mq5); AIC=AIC[,-1]
BIC=BIC(mq1,mq2,mq3,mq4, mq5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1q,s2q,s3q,s4q,s5q, a1q,a2q,a3q,a4q, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)
    }
# linear plateau
    f3=function(data){
                      ml = lm(data[, 3] ~ data[, 1])
                      mq = lm(data[, 3] ~ data[, 1] + I(data[, 1]^2))
                      c1 = coef(ml)[[1]]
                      c2 = coef(ml)[[2]]
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
                      pc = -0.5 * c4/c5
			ff=function(x){c3+c4*x+c5*x^2}
			pp=ff(pc)
			a=start[1]
			b=start[2]
			c=start[3]
			a11=ifelse(a==1,pp,a)
			b11=ifelse(b==1,c2,b)
			c11=ifelse(c==1,pc,c)
		        
# normal
m1 <- nls(y ~ a + b * (x - c) * (x <= c), start = list(a = a11, 
                                                                        b = b11, c = c11), data = data, control = nls.control(maxiter = 6000))
# a test
m2 <- nls(y ~ a[treatment] + b * (x - c) * (x <= c), start = list(a = c(a11,a11),b = b11, c = c11), data = data, control = nls.control(maxiter = 6000))
# b test
m3 <- nls(y ~ a + b[treatment] * (x - c) * (x <= c), start = list(a = c(a11),b = c(b11,b11), c = c11), data = data, control = nls.control(maxiter = 6000))
# c test
m4 <- nls(y ~ a + b * (x - c[treatment]) * (x <= c[treatment]), start = list(a = c(a11),b = c(b11), c = c(c11,c11)), data = data, control = nls.control(maxiter = 6000))
# model test
m5 <- nls(y ~ a[treatment] + b[treatment] * (x - c[treatment]) * (x <= c[treatment]), start = list(a = c(a11,a11),b = c(b11,b11), c = c(c11,c11)), data = data, control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)
    }  
    # quadratic plateau
    f4=function(data){
                      mq = lm(data[, 3] ~ data[, 1] + I(data[, 1]^2))
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
			a=start[1]
			b=start[2]
			c=start[3]
			a11=ifelse(a==1,c3,a)
			b11=ifelse(b==1,c4,b)
			c11=ifelse(c==1,c5,c)

# normal
m1 <- nls(y ~ (a + b * x + c * I(x^2)) * (x <= -0.5 *b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c),  start = list(a = a11, b = b11, c = c11), data = data,control = nls.control(maxiter = 6000))
# a test
m2 <- nls(y ~ (a[treatment] + b * x + c * I(x^2)) * (x <= -0.5 *b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c),  start = list(a = c(a11,a11), b = b11, c = c11), data = data,control = nls.control(maxiter = 6000))
# b test
m3 <- nls(y ~ (a + b[treatment] * x + c * I(x^2)) * (x <= -0.5 *b[treatment]/c) + (a + I(-b[treatment]^2/(4 * c))) * (x > -0.5 * b[treatment]/c),  start = list(a = c(a11), b = c(b11,b11), c = c11), data = data,control = nls.control(maxiter = 6000))
# c test
m4 <- nls(y ~ (a + b * x + c[treatment] * I(x^2)) * (x <= -0.5 *b/c[treatment]) + (a + I(-b^2/(4 * c[treatment]))) * (x > -0.5 * b/c[treatment]),  start = list(a = c(a11), b = c(b11), c = c(c11,c11)), data = data,control = nls.control(maxiter = 6000))
# model test
m5 <- nls(y ~ (a[treatment] + b[treatment] * x + c[treatment] * I(x^2)) * (x <= -0.5 *b[treatment]/c[treatment]) + (a[treatment] + I(-b[treatment]^2/(4 * c[treatment]))) * (x > -0.5 * b[treatment]/c[treatment]),  start = list(a = c(a11,a11), b = c(b11,b11), c = c(c11,c11)), data = data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)
    } 
    # bi linear
    f5=function(data){   
s=start                   
fp=function(a,b,c,x,d){ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)}                   
# normal
m1=nls(y~fp(a,b,c,x,d), start=list(a=s[1],b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# a test
m2=nls(y~fp(a[treatment],b,c,x,d), start=list(a=c(s[1],s[1]),b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~fp(a,b[treatment],c,x,d), start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# c test
m4=nls(y~fp(a,b,c[treatment],x,d), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3]), d=s[4]), data=data,control = nls.control(maxiter = 6000))
# d test
m5=nls(y~fp(a,b,c,x,d[treatment]), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
# model test
m6=nls(y~fp(a[treatment],b[treatment],c[treatment],x,d[treatment]), start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3], s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
s6=summary(m6)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
a5=anova(m1,m6)
AIC=AIC(m1,m2,m3,m4, m5,m6); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5,m6); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2")
l2=list(s1,s2,s3,s4,s5,s6, a1,a2,a3,a4,a5, par)
names(l2)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2","test a parameter", "test b parameter","test c parameter","test d parameter","test model", "AIC and BIC") 
return(l2)
    }
    # exponential
    f6=function(data){
s=start
# normal
m1=nls(y~a*exp(b*x) ,start=list(a=s[1],b=s[2]),data=data,control = nls.control(maxiter = 6000)) 
# a test
m2=nls(y~a[treatment]*exp(b*x) ,start=list(a=c(s[1],s[1]),b=s[2]),data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*exp(b[treatment]*x) ,start=list(a=c(s[1]),b=c(s[2],s[2])),data=data,control = nls.control(maxiter = 6000))
# model test
m4=nls(y~a[treatment]*exp(b[treatment]*x) ,start=list(a=c(s[1],s[1]),b=c(s[2],s[2])),data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
AIC=AIC(m1,m2,m3,m4); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2")
l1=list(s1,s2,s3,s4,a1,a2,a3, par)
names(l1)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2","test a parameter", "test b parameter","test model", "AIC and BIC") 
return(l1)
    } 
    
# logistic model
 f7=function(data){
s=start
# normal
m1=nls(y~a*(1+b*(exp(-c*x)))^-1, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# a test
m2=nls(y~a[treatment]*(1+b*(exp(-c*x)))^-1, data=data, start=list(a=c(s[1],s[1]),b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*(1+b[treatment]*(exp(-c*x)))^-1, data=data, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3]),control = nls.control(maxiter = 6000))
# c test
m4=nls(y~a*(1+b*(exp(-c[treatment]*x)))^-1, data=data, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
# model test
m5=nls(y~a[treatment]*(1+b[treatment]*(exp(-c[treatment]*x)))^-1, data=data, start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)
    }
    
# Von Bertalanffy model
f8=function(data){
s=start
# normal
m1=nls(y~a*(1-b*(exp(-c*x)))^3, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# a test
m2=nls(y~a[treatment]*(1-b*(exp(-c*x)))^3, data=data, start=list(a=c(s[1],s[1]),b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*(1-b[treatment]*(exp(-c*x)))^3, data=data, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3]),control = nls.control(maxiter = 6000))
# c test
m4=nls(y~a*(1-b*(exp(-c[treatment]*x)))^3, data=data, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
# model test
m5=nls(y~a[treatment]*(1-b[treatment]*(exp(-c[treatment]*x)))^3, data=data, start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)
    }  
# Brody model
f9=function(data){
s=start
# normal
m1=nls(y~a*(1-b*(exp(-c*x))), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
 # a test
m2=nls(y~a[treatment]*(1-b*(exp(-c*x))), data=data, start=list(a=c(s[1],s[1]),b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*(1-b[treatment]*(exp(-c*x))), data=data, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3]),control = nls.control(maxiter = 6000))
# c test
m4=nls(y~a*(1-b*(exp(-c[treatment]*x))), data=data, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
# model test
m5=nls(y~a[treatment]*(1-b[treatment]*(exp(-c[treatment]*x))), data=data, start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)                 
    }
# Gompertz model
f10=function(data){
s=start
# normal
m1=nls(y~a*exp(-b*exp(-c*x)), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000)) 
# a test
m2=nls(y~a[treatment]*exp(-b*exp(-c*x)), data=data, start=list(a=c(s[1],s[1]),b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*exp(-b[treatment]*exp(-c*x)), data=data, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3]),control = nls.control(maxiter = 6000))
# c test
m4=nls(y~a*exp(-b*exp(-c[treatment]*x)), data=data, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
# model test
m5=nls(y~a[treatment]*exp(-b[treatment]*exp(-c[treatment]*x)), data=data, start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)                    
    }
    
    # lactation curve
    f11=function(data){
s=start
# normal
m1=nls(y~(a*x^b)*exp(-c*x), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# a test
m2=nls(y~(a[treatment]*x^b)*exp(-c*x), data=data, start=list(a=c(s[1],s[1]),b=s[2],c=s[3]),control = nls.control(maxiter = 6000))
# b test
m3=nls(y~(a*x^b[treatment])*exp(-c*x), data=data, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3]),control = nls.control(maxiter = 6000))
# c test
m4=nls(y~(a*x^b)*exp(-c[treatment]*x), data=data, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
# model test
m5=nls(y~(a[treatment]*x^b[treatment])*exp(-c[treatment]*x), data=data, start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3])),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)        
                  
                       # a=17;b=0.25; c=0.004
    }
    
    y1=c(25,24,26,28,30,31,27,26,25,24,23,24,22,21,22,20,21,19,18,17,18,18,16,17,15,16,14)
    y2=y1+2
    y=c(y1,y2)
    x1=c(15,15,15,75,75,75,135,135,135,195,195,195,255,255,255,315,315,315,375,375,375,435,435,435,
         495,495,495);x2=x1
    x=c(x1,x2); treatment=rep(c("t1","t2"), each=27)
    
    dados=data.frame(x,treatment,y)
    
    
    # ruminal degradation curve
    f12=function(data){
s=start
# normal
m1=nls(y ~ a + b * (1 - exp(-c * x)), data=data, start=list(a=20,b=60,c=0.05),control = nls.control(maxiter = 6000))
# a test
m2=nls(y ~ a[treatment] + b * (1 - exp(-c * x)), data=data, start=list(a=c(20,20),b=60,c=0.05),control = nls.control(maxiter = 6000))
# b test
m3=nls(y ~ a + b[treatment] * (1 - exp(-c * x)), data=data, start=list(a=c(20),b=c(60,60),c=0.05),control = nls.control(maxiter = 6000))
# c test
m4=nls(y ~ a + b* (1 - exp(-c[treatment]  * x)), data=data, start=list(a=c(20),b=c(60),c=c(0.05,0.05)),control = nls.control(maxiter = 6000))
# d test
m5=nls(y ~ a[treatment] + b[treatment]* (1 - exp(-c[treatment]  * x)), data=data, start=list(a=c(20,20),b=c(60,60),c=c(0.05,0.05)),control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
AIC=AIC(m1,m2,m3,m4, m5); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2")
l2=list(s1,s2,s3,s4,s5, a1,a2,a3,a4, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test model", "AIC and BIC") 
return(l2)        

    }
    
    tempo=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
    gas=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    
    d=data.frame(tempo,gas)
    
    m=nls(gas~(vf1/(1+exp(2-4*k1*(tempo-l))))+(vf2/(1+exp(2-4*k2*(tempo-l)))), start=list(vf1=19,vf2=4,k1=0.05,k2=0.005,l=5), data=d)


    t1=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
    t2=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
    tempo=c(t1,t2)
    gas1=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    gas2=gas1*0.95
    gas=c(gas1,gas2)
    treatment=rep(c("t1","t2"), each=14)
    
    d=data.frame(tempo,treatment,gas)
    
    
    # logistico bicompartimental
    
    f13=function(data){
s=start
# normal
m1=nls(y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e)))), start=list(a=s[1],b=s[2],c=s[3],d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000))
# a test
m2=nls(y~(a[treatment]/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e)))), start=list(a=c(s[1],s[1]),b=s[2],c=s[3],d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~(a/(1+exp(2-4*c*(x-e))))+(b[treatment]/(1+exp(2-4*d*(x-e)))), start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3],d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000))
# c test
m4=nls(y~(a/(1+exp(2-4*c[treatment]*(x-e))))+(b/(1+exp(2-4*d*(x-e)))), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3]),d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000))
# d test
m5=nls(y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d[treatment]*(x-e)))), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3]),d=c(s[4],s[4]),e=s[5]), data=data,control = nls.control(maxiter = 6000))
# e test
m6=nls(y~(a/(1+exp(2-4*c*(x-e[treatment]))))+(b/(1+exp(2-4*d*(x-e[treatment])))), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3]),d=c(s[4]),e=c(s[5],s[5])), data=data,control = nls.control(maxiter = 6000))
# model test
m7=nls(y~(a[treatment]/(1+exp(2-4*c[treatment]*(x-e[treatment]))))+(b[treatment]/(1+exp(2-4*d[treatment]*(x-e[treatment])))), start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3],s[3]),d=c(s[4],s[4]),e=c(s[5],s[5])), data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
s6=summary(m6)[[10]]
s7=summary(m7)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
a5=anova(m1,m6)
a6=anova(m1,m7)
AIC=AIC(m1,m2,m3,m4, m5, m6, m7); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5, m6, m7); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c d e","a1 a2 b c d e","a b1 b2 c d e","a b c1 c2 d e","a b c d1 d2 e","a b c d e1 e2","a1 a2 b1 b2 c1 c2 d1 d2 e1 e2")
l2=list(s1,s2,s3,s4,s5,s6,s7, a1,a2,a3,a4,a5,a6, par)
names(l2)=c("a b c d e","a1 a2 b c d e","a b1 b2 c d e","a b c1 c2 d e","a b c d1 d2 e","a b c d e1 e2","a1 a2 b1 b2 c1 c2 d1 d2 e1 e2","test a parameter", "test b parameter","test c parameter","test d parameter","test e parameter","test model", "AIC and BIC") 
return(l2)




    }


# exponential (allometric model)
    f14=function(data){
s=start
# normal
m1=nls(y~a*(x^b) ,start=list(a=s[1],b=s[2]),data=data,control = nls.control(maxiter = 6000)) 
# a test
m2=nls(y~a[treatment]*(x^b) ,start=list(a=c(s[1],s[1]),b=s[2]),data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a*(x^b[treatment]) ,start=list(a=c(s[1]),b=c(s[2],s[2])),data=data,control = nls.control(maxiter = 6000))
# model test
m4=nls(y~a[treatment]*(x^b[treatment]) ,start=list(a=c(s[1],s[1]),b=c(s[2],s[2])),data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
AIC=AIC(m1,m2,m3,m4); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2")
l1=list(s1,s2,s3,s4,a1,a2,a3, par)
names(l1)=c("a b","a1 a2 b","a b1 b2","a1 a2 b1 b2","test a parameter", "test b parameter","test model", "AIC and BIC") 
return(l1)
    } 

# cubic
f15=function(data){     
mq1=nls(y~a+b*x+c*x^2+d*x^3, data=data, start=c(a=1,b=1,c=1, d=1))
# a
mq2=nls(y~a[treatment]+b*x+c*x^2+d*x^3, data=data, start=list(a=c(1,1),b=c(1), c=c(1), d=c(1)))
# b
mq3=nls(y~a+b[treatment]*x+c*x^2+d*x^3, data=data, start=list(a=c(1),b=c(1,1), c=c(1), d=c(1)))
# c
mq4=nls(y~a+b*x+c[treatment]*x^2+d*x^3, data=data, start=list(a=c(1),b=c(1), c=c(1,1), d=c(1)))
# d
mq5=nls(y~a+b*x+c*x^2+d[treatment]*x^3, data=data, start=list(a=c(1),b=c(1), c=c(1), d=c(1,1)))
# model
mq6=nls(y~a[treatment]+b[treatment]*x+c[treatment]*x^2+d[treatment]*x^3, data=data, start=list(a=c(1,1),b=c(1,1), c=c(1,1), d=c(1,1)))
s1q=summary(mq1)[[10]]
s2q=summary(mq2)[[10]]
s3q=summary(mq3)[[10]]
s4q=summary(mq4)[[10]]
s5q=summary(mq5)[[10]]
s6q=summary(mq5)[[10]]
a1q=anova(mq1,mq2)
a2q=anova(mq1,mq3)
a3q=anova(mq1,mq4)
a4q=anova(mq1,mq5)
a5q=anova(mq1,mq6)
AIC=AIC(mq1,mq2,mq3,mq4, mq5, mq6); AIC=AIC[,-1]
BIC=BIC(mq1,mq2,mq3,mq4, mq5, mq6); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2")
l2=list(s1q,s2q,s3q,s4q,s5q,s6q, a1q,a2q,a3q,a4q,a5q, par)
names(l2)=c("a b c","a1 a2 b c","a b1 b2 c","a b c1 c2","a b c d1 d2","a1 a2 b1 b2 c1 c2","test a parameter", "test b parameter","test c parameter","test d parameter","test model", "AIC and BIC") 
return(l2)
    }
    

# richards model
    f16=function(data){   
s=start                   
# normal
m1=nls(y~a/(1+b*(exp(-c*x)))^d, start=list(a=s[1],b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# a test
m2=nls(y~a[treatment]/(1+b*(exp(-c*x)))^d, start=list(a=c(s[1],s[1]),b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~a/(1+b[treatment]*(exp(-c*x)))^d, start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# c test
m4=nls(y~a/(1+b*(exp(-c[treatment]*x)))^d, start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3]), d=s[4]), data=data,control = nls.control(maxiter = 6000))
# d test
m5=nls(y~a/(1+b*(exp(-c*x)))^d[treatment], start=list(a=c(s[1]),b=c(s[2]),c=c(s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
# model test
m6=nls(y~a[treatment]/(1+b[treatment]*(exp(-c[treatment]*x)))^d[treatment], start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3], s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
s6=summary(m6)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
a5=anova(m1,m6)
AIC=AIC(m1,m2,m3,m4, m5,m6); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5,m6); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2")
l2=list(s1,s2,s3,s4,s5,s6, a1,a2,a3,a4,a5, par)
names(l2)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2","test a parameter", "test b parameter","test c parameter","test d parameter","test model", "AIC and BIC") 
return(l2)
    }

# schnute model
    f17=function(data){   
s=start  
t1=min(data[,2],na.rm = TRUE)
t2=max(data[,2],na.rm = TRUE)                 
# normal
m1=nls(y~(a^d+ ((b^d)-(a^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d), start=list(a=s[1],b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# a test
m2=nls(y~(a[treatment]^d+ ((b^d)-(a[treatment]^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d), start=list(a=c(s[1],s[1]),b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# b test
m3=nls(y~(a^d+ ((b[treatment]^d)-(a^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d), start=list(a=c(s[1]),b=c(s[2],s[2]),c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
# c test
m4=nls(y~(a^d+ ((b^d)-(a^d) )*((1-exp(-c[treatment]*(x-t1)))/ (1-exp(-c[treatment]*(t2-t1)))))^(1/d), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3],s[3]), d=s[4]), data=data,control = nls.control(maxiter = 6000))
# d test
m5=nls(y~(a^d[treatment]+ ((b^d[treatment])-(a^d[treatment]) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d[treatment]), start=list(a=c(s[1]),b=c(s[2]),c=c(s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
# model test
m6=nls(y~(a[treatment]^d[treatment]+ ((b[treatment]^d[treatment])-(a[treatment]^d[treatment]) )*((1-exp(-c[treatment]*(x-t1)))/ (1-exp(-c[treatment]*(t2-t1)))))^(1/d[treatment]), start=list(a=c(s[1],s[1]),b=c(s[2],s[2]),c=c(s[3], s[3]), d=c(s[4],s[4])), data=data,control = nls.control(maxiter = 6000))
s1=summary(m1)[[10]]
s2=summary(m2)[[10]]
s3=summary(m3)[[10]]
s4=summary(m4)[[10]]
s5=summary(m5)[[10]]
s6=summary(m6)[[10]]
a1=anova(m1,m2)
a2=anova(m1,m3)
a3=anova(m1,m4)
a4=anova(m1,m5)
a5=anova(m1,m6)
AIC=AIC(m1,m2,m3,m4, m5,m6); AIC=AIC[,-1]
BIC=BIC(m1,m2,m3,m4, m5,m6); BIC=BIC[,-1]
par=data.frame(AIC,BIC)
rownames(par)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2")
l2=list(s1,s2,s3,s4,s5,s6, a1,a2,a3,a4,a5, par)
names(l2)=c("a b c d","a1 a2 b c d","a b1 b2 c d","a b c1 c2 d","a b c d1 d2","a1 a2 b1 b2 c1 c2 d1 d2","test a parameter", "test b parameter","test c parameter","test d parameter","test model", "AIC and BIC") 
return(l2)
    }

 
    fun=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17)
    
    fu=fun[[model]]
    r1=lapply(r2, fu)

    name=levels(data$treatment)
    c=combn(name,2)

	io=1:ncol(c)
	fc=function(io){
	p=paste(c[,io], collapse=" vs ")
	return(p)
	}

lc=lapply(io,fc)
comb=unlist(lc)

names(r1)=comb
  
    return(r1)
    
}
