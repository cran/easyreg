	er1 <-
function(data, model=1, start=c(a=1,b=1,c=1,d=1,e=1), mixed=FALSE, digits=6, alpha=0.05){
    s=start
    mixed=mixed

tmodel<-function(m, data){
i=if(class(m)=="lm") 1 else 2
ll=list(lm(data[,2]~1),lm(y~1, data=data))
m0=ll[[i]]
a=anova(m,m0)
df=c(abs(a[,3][2]),a[,1][1],a[,1][2])
sq=c(abs(a[,4][2]),a[,2][1],a[,2][2])
qm=sq/df
d1=data.frame(df,sq,qm);d1=round(d1,4);
f=qm[1]/qm[2]
fc=c(round(f,4),"-","-")
p=pf(f,df[1],df[2], lower.tail=FALSE)
p=if(p<0.001) "<0.001" else p=round(p,4)
ps=c(p,"-","-")
source=c("regression", "residuals", "total")
res=data.frame(d1,fc,ps)
colnames(res)=c("df","sum of squares", "mean squares", "Fcal", "p-value")
rownames(res)=source
return(res)
}

fwm<-function(m){
s=as.numeric(summary(m)[6])
v=s^2
a=anova(m)[-1,]
gle=a[,2][1]
glr=nrow(a)
sqr=sum(a[,3]*v)
qmr=sqr/glr
fcr=qmr/v
p=pf(fcr,glr, gle, lower.tail = FALSE)
r=data.frame(Df=round(c(glr,gle),4), SS=round(c(sqr, v*gle),4), Mean_Square=round(c(qmr,v),4), Fcal=c(round(fcr,2), "-"),p_value=c(round(p,4), "-"))
rownames(r)=c("Whole model", "Residuals")
return(r)
}

fcon<-function(m){
s=diag(vcov(m))^0.5
gl=list(summary(m)$df[2],m$df.residual, m$fixDF[[1]][[1]])
i=if(class(m)[1]=="nls") 1
i=if(class(m)[1]=="lm") 2 else i
i=if(class(m)[1]=="nlme") 3 else i
i=if(class(m)[1]=="lme") 3 else i
gl=gl[[i]]
t=abs(qt((alpha/2),gl))
f1=function(x)fixef(x);f2=function(x) coef(x)
cc=list(f1,f2)
i=if(class(m)[1]=="nls") 2
i=if(class(m)[1]=="lm") 2 else i
i=if(class(m)[1]=="nlme") 1 else i
i=if(class(m)[1]=="lme") 1 else i
c=as.numeric(cc[[i]](m))
di=t*s
c1=c-di
c2=c+di
l=letters[1:length(coef(m))]
r=data.frame(l,c,s,c1,c2)	
p1=paste((alpha/2)*100,"%",sep="")
p2=paste((1-(alpha/2))*100,"%",sep="")
p1=paste("IC",p1, sep=" ")
p2=paste("IC",p2, sep=" ")
r=data.frame(Parameters=r[,1],round(r[,-1], digits))
names(r)=c("parameter", "estimate", "standard error", p1,p2) 
rownames(r)=NULL
return(r)
}

  opt1=function(data){
    d1 = as.data.frame(data[,1])
    d2 = as.data.frame(data[,-1])
    f = function(h) {
        data.frame(d1, d2[h])}
    h = length(d2)
    h = 1:h
    l = lapply(h, f)
    return(l)
    }

    opt2=function(data){
    d3 = as.data.frame(data[,c(1,2)])
    d4 = as.data.frame(data[,-c(1,2)])
    ff = function(h) {
        data.frame(d3, d4[h])}
    hh = length(d4)
    hh = 1:hh
    ll = lapply(hh, ff)
    return(ll)
    }


 

cv <- function(x) {
        sd = (deviance(x)/df.residual(x))^0.5
        mm = mean(fitted(x))
        r = 100 * sd/mm
        return(r)
    }

cvmm <- function(x) {
	sd = x$sigma
        mm = mean(fitted(x))
        r = 100 * sd/mm
        return(r)
    }
        

	fr=function(m){
	mixed=mixed
        r=resid(m)
	t=1:length(r)
	i=ifelse(length(r)>5000, 2,1)
	jr=function(r,aa)r+aa-aa
	jsample=function(r,aa)sample(r,aa)
	rr=list(jr,jsample)
	rr=rr[[i]](r,5000)
        s <- shapiro.test(rr)
	i=ifelse(mixed==FALSE,1,2)
	lf=list(cv,cvmm)
	FUNN=lf[[i]]
        cvf = FUNN(m)
	names(r)=t
        rd=as.data.frame((sort(sqrt(r^2),decreasing=TRUE)))
        rl=as.list(rownames(rd))
        r1=rl[[1]];r2=rl[[2]];r3=rl[[3]]
        d=c(s$"p.value",cvf,as.numeric(r1),as.numeric(r2),as.numeric(r3))
        return(d)}

pres=function(m,name="name"){
	name=name
	r=resid(m)
	r=scale(r)
	t=1:length(r)
	ll=labels=1:length(r);mo=median(ll)
	g1=function(r){plot(r~t, pch="",ylim=c(-4,4), ylab="Standardized residuals", xlab="Sequence data", 		main="Standardized residuals vs Sequence data",axes=FALSE);axis(2,c(-4,-3.5,-3,-2.5,-2,-1,0,1,2,2.5,3,3.5,4));abline(h=2.5, lty=2);abline(h=-2.5,lty=2);abline(h=3.5, lty=2, col=2);abline(h=-3.5,lty=2, col=2); text(2.5,2.7, "2.5 z-score");text(2.5,-2.7, "-2.5 z-score");text(2.5,3.7, "3.5 z-score");text(2.5,-3.7, "-3.5 z-score");text(t,r,labels=1:length(r));text(mo,4,name, cex=2, col=2)}
	g1(r)
	}



R2 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        return(r1)
    }

R3 <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
        r1 <- (sqt - deviance(m))/sqt
        p1 <- (gl/((gl + 1) - (length(coef(m) + 1))) * (1 - 
            r1))
        r2 <- 1 - p1; return(r2)
    }

R2mm <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
	dev=sum(resid(m)^2)
        r1 <- (sqt - dev)/sqt
        return(r1)
    }


R3mm <- function(m) {
        gl <- length(fitted(m)) - 1
        sqt <- var((fitted(m) + resid(m))) * gl
	dev=sum(resid(m)^2)
        r1 <- (sqt - dev)/sqt
        p1 <- (gl/((gl + 1) - (length(fixef(m) + 1))) * (1 - 
            r1))
        r2 <- 1 - p1; return(r2)
    }
    
    # linear
    f1=function(data){names(data) = c("x", "y") 
                      ml = lm(data[,2] ~ data[, 1])
                      c1 = coef(ml)[[1]]
                      c2 = coef(ml)[[2]]
                      c=coef(ml); s=summary(ml); a=c[1];b=c[2]; fr1=fr(ml); 
re1=resid(ml);re2=scale(re1)
con=fcon(ml)
l=c(a,b,summary(ml)$coefficients[,4],(summary(ml)$sigma)^2, R2(ml), R3(ml),AIC(ml), BIC(ml),fr1) 
l=round(l,digits);l=as.data.frame(l); names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "p-value t.test for a", "p-value t.test for b","residual mean square", 
                                    "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(ml, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Linear model","Parameters, standard error and confidence intervals","Whole model test","Residuals of linear model", "Residuals standardized of linear model") 
return(l1)
                      }
    
    
    # quadratic
    f2=function(data){names(data) = c("x", "y") 
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
                      c=coef(mq); s=summary(mq); a=c[1];b=c[2];c=c[3];pm = a - (b^2)/(4 * c)
                      pc = -0.5 * b/c
fr1=fr(mq)
re1=resid(mq);re2=scale(re1)
con=fcon(mq)
 l=c(a,b,c,summary(mq)$coefficients[,4], (summary(mq)$sigma)^2, R2(mq), R3(mq),AIC(mq), BIC(mq), pm,pc,fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b",
                                    "p-value t.test for c","residual mean square", "r-squared", "adjusted r-squared","AIC", "BIC","maximum or minimum value for y", "critical point in x", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(mq, data)
                    l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Quadratic model","Parameters, standard error and confidence intervals","Whole model test","Residuals of quadratic model", "Residuals standardized of quadratic model") 
return(l1)
    }
    
   # linear plateau
    f3=function(data){names(data) = c("x", "y")
                      ml = lm(data[, 2] ~ data[, 1])
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
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
		        m <- nls(y ~ a + b * (x - c) * (x <= c), start = list(a = a11, 
                                                                        b = b11, c = c11), data = data, control = nls.control(maxiter = 6000))
                      c=coef(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
                      l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), a,c, fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "residual mean square","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","maximum or minimum value for y", "critical point in x", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(m, data)
                       l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Linear plateau model","Parameters, standard error and confidence intervals","Whole model test","Residuals of linear plateau model", "Residuals standardized of linear plateau model") 
return(l1)
    }
    
    
    
    # quadratic plateau
    f4=function(data){names(data) = c("x", "y")
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2))
                      c3 = coef(mq)[[1]]
                      c4 = coef(mq)[[2]]
                      c5 = coef(mq)[[3]]
			
			a=start[1]
			b=start[2]
			c=start[3]
			a11=ifelse(a==1,c3,a)
			b11=ifelse(b==1,c4,b)
			c11=ifelse(c==1,c5,c)


                      m <- nls(y ~ (a + b * x + c * I(x^2)) * (x <= -0.5 * 
                                                                   b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c), 
                               start = list(a = a11, b = b11, c = c11), data = data, 
                               control = nls.control(maxiter = 6000))
                      c=coef(m); a=c[1];b=c[2];c=c[3]
                      pmm = (a + I(-b^2/(4 * c)))
                      pcc = -0.5 * b/c
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
                      l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), pmm,pcc, fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c","residual mean square", "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","maximum or minimum value for y", "critical point in x", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(m, data) 
                     l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Quadratic plateau model","Parameters, standard error and confidence intervals","Whole model test","Residuals of quadratic plateau model", "Residuals standardized of quadratic plateau model") 
return(l1)
    }
    
    # bi linear
    f5=function(data){names(data) = c("x", "y")     
                      fp=function(a,b,c,x,d){ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)}
                      m=nls(y~fp(a,b,c,x,d), start=list(a=s[1],b=s[2],c=s[3], d=s[4]), data=data,control = nls.control(maxiter = 6000))
c=coef(m); a=c[1];b=c[2];c1=c;c=c[3];d=c1[4]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
                      l=c(a,b,c,d,summary(m)[11][[1]][13], summary(m)[11][[1]][14], summary(m)[11][[1]][15],summary(m)[11][[1]][16], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c","coefficient d", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "p-value t.test for d","residual mean square","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Broken line model","Parameters, standard error and confidence intervals","Whole model test","Residuals of broken line model", "Residuals standardized of broken line model") 
return(l1)

    }
    
    # exponential
    f6=function(data){names(data) = c("x", "y")
                      m=nls(y~a*exp(b*x) ,start=list(a=s[1],b=s[2]),data=data,control = nls.control(maxiter = 6000));c=coef(m); a=c[1];b=c[2]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
                      l=c(a,b,summary(m)[11][[1]][7], summary(m)[11][[1]][8],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "p-value t.test for a", "p-value t.test for b", "residual mean square",
                                    "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Exponential model","Parameters, standard error and confidence intervals","Whole model test","Residuals of exponential model", "Residuals standardized of exponential model") 
return(l1)
    } 
    
    # logistic model
    f7=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1+b*(exp(-c*x)))^-1, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "residual mean square","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Logistic model","Parameters, standard error and confidence intervals","Whole model test","Residuals of logistic model", "Residuals standardized of logistic model") 
return(l1)
    }
    
    # Von Bertalanffy model
    f8=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1-b*(exp(-c*x)))^3, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
 l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "residual mean square","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(m, data) 
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Von Bertalanffy model","Parameters, standard error and confidence intervals","Whole model test","Residuals of Von Bertalanffy model", "Residuals standardized of Von Bertalanffy model") 
return(l1)
    }
    
    # Brody model
    f9=function(data){names(data) = c("x", "y") 
                      m=nls(y~a*(1-b*(exp(-c*x))), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                    "p-value t.test for c", "residual mean square","r-squared", "adjusted r-squared", 
                                    "AIC", "BIC","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Brody model","Parameters, standard error and confidence intervals","Whole model test","Residuals of Brody model", "Residuals standardized of Brody model") 
return(l1)
    }
    
    # Gompertz model
    f10=function(data){names(data) = c("x", "y") 
                       m=nls(y~a*exp(-b*exp(-c*x)), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c", "residual mean square","r-squared", "adjusted r-squared", 
                                     "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Gompertz model","Parameters, standard error and confidence intervals","Whole model test","Residuals of Gompertz model", "Residuals standardized of Gompertz model") 
return(l1)
    }
    
    # lactation curve
    f11=function(data){names(data) = c("x", "y") 
                       m=nls(y~(a*x^b)*exp(-c*x), data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]; pm=a*(b/c)^b*exp(-b);pc=b/c
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m), pm,pc, fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c","residual mean square", "r-squared", "adjusted r-squared", 
                                     "AIC", "BIC","maximum or minimum value for y", "critical point in x","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Lactation model","Parameters, standard error and confidence intervals","Whole model test","Residuals of lactation model", "Residuals standardized of lactation model") 
return(l1)
                       # a=17;b=0.25; c=0.004
    }
    
    y1=c(25,24,26,28,30,31,27,26,25,24,23,24,22,21,22,20,21,19,18,17,18,18,16,17,15,16,14)
    x1=c(15,15,15,75,75,75,135,135,135,195,195,195,255,255,255,315,315,315,375,375,375,435,435,435,
         495,495,495)
    
    dados=data.frame(x1,y1)
    
    
    # ruminal degradation curve
    f12=function(data){names(data) = c("x", "y") 
                       m=nls(y ~ a + b * (1 - exp(-c * x)), data=data, start=list(a=20,b=60,c=0.05),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
 l=c(a,b,c,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m),fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                       rownames(l)=c("coefficient a", "coefficient b", 
                                     "coefficient c", "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c","residual mean square", "r-squared", "adjusted r-squared", 
                                     "AIC", "BIC","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Ruminal degradation model","Parameters, standard error and confidence intervals","Whole model test","Residuals of ruminal degradation model", "Residuals standardized of ruminal degradation model") 
return(l1)
    }
    
    tempo=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
    gas=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    d=data.frame(tempo,gas)
    m=nls(gas~(vf1/(1+exp(2-4*k1*(tempo-l))))+(vf2/(1+exp(2-4*k2*(tempo-l)))), start=list(vf1=19,vf2=4,k1=0.05,k2=0.005,l=5), data=d)
    
    # logistico bicompartimental   
    f13=function(data){names(data) = c("x", "y") 
                       m=nls(y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e)))), start=list(a=s[1],b=s[2],c=s[3],d=s[4],e=s[5]), data=data,control = nls.control(maxiter = 6000));cc=coef(m); s=summary(m); a=cc[1];b=cc[2];c=cc[3];d=cc[4];e=cc[5]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
 l=c(a,b,c,d,e,summary(m)[11][[1]][16], summary(m)[11][[1]][17], summary(m)[11][[1]][18],summary(m)[11][[1]][19],summary(m)[11][[1]][20], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m),fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                       rownames(l)=c("coefficient a", "coefficient b", "coefficient c","coefficient d","coefficient e",
                                     "p-value t.test for a", "p-value t.test for b", 
                                     "p-value t.test for c","p-value t.test for d","p-value t.test for e" ,"residual mean square","r-squared", "adjusted r-squared", 
                                     "AIC", "BIC","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Logistic bi-compartmental model","Parameters, standard error and confidence intervals","Whole model test","Residuals of logistic bi-compartmental model", "Residuals standardized of logistic bi-compartmental model") 
return(l1)
    }

# exponential (allometric)
    f14=function(data){names(data) = c("x", "y")
                      m=nls(y~a*(x^b) ,start=list(a=s[1],b=s[2]),data=data,control = nls.control(maxiter = 6000));c=coef(m); a=c[1];b=c[2]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
                      l=c(a,b,summary(m)[11][[1]][7], summary(m)[11][[1]][8], (summary(m)$sigma)^2,R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "p-value t.test for a", "p-value t.test for b","residual mean square", 
                                    "r-squared", "adjusted r-squared", 
                                    "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Exponential (allometric) model","Parameters, standard error and confidence intervals","Whole model test","Residuals of allometric model", "Residuals standardized of allometric model") 
return(l1)
    } 

# cubic
    f15=function(data){names(data) = c("x", "y") 
                      mq = lm(data[, 2] ~ data[, 1] + I(data[, 1]^2)+I(data[, 1]^3))
                      c=coef(mq); s=summary(mq); a=c[1];b=c[2];d=c[4];c=c[3]
fr1=fr(mq)
re1=resid(mq);re2=scale(re1)
con=fcon(mq)
 l=c(a,b,c,d,summary(mq)$coefficients[,4], (summary(mq)$sigma)^2, R2(mq), R3(mq),AIC(mq), BIC(mq),fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c","coefficient d", "p-value t.test for a", "p-value t.test for b",
                                    "p-value t.test for c","p-value t.test for d","residual mean square", "r-squared", "adjusted r-squared","AIC", "BIC","p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(mq, data)
                    l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Cubic model","Parameters, standard error and confidence intervals","Whole model test","Residuals of cubic model", "Residuals standardized of cubic model") 
return(l1)
    }
    
 # richards model
    f16=function(data){names(data) = c("x", "y") 
                      m=nls(y~a/((1+b*exp( -c*x) )^d), data=data, start=list(a=s[1],b=s[2],c=s[3], d=s[4]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];d=c[4];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,d,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],summary(m)[11][[1]][13],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c","coefficient d", "p-value t.test for a", "p-value t.test for b", "p-value t.test for c", "p-value t.test for d","residual mean square","r-squared", "adjusted r-squared","AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Richards model","Parameters, standard error and confidence intervals","Whole model test","Residuals of richards model", "Residuals standardized of richards model") 
return(l1)
    }

# schnute model
    f17=function(data){names(data) = c("x", "y") 
t1=min(data$y, na.rm = TRUE)
t2=max(data$y, na.rm = TRUE)
                      m=nls(y~(a^d+ ((b^d)-(a^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d),  data=data, start=list(a=s[1],b=s[2],c=s[3], d=s[4]),control = nls.control(maxiter = 6000));c=coef(m); s=summary(m); a=c[1];b=c[2];d=c[4];c=c[3]
fr1=fr(m)
con=fcon(m)
re1=resid(m);re2=scale(re1)
l=c(a,b,c,d,summary(m)[11][[1]][10], summary(m)[11][[1]][11], summary(m)[11][[1]][12],summary(m)[11][[1]][13],(summary(m)$sigma)^2, R2(m), R3(m),AIC(m), BIC(m), fr1); l=round(l,digits);l=as.data.frame(l)
names(l)="estimates"
                      rownames(l)=c("coefficient a", "coefficient b", 
                                    "coefficient c","coefficient d", "p-value t.test for a", "p-value t.test for b", "p-value t.test for c", "p-value t.test for d","residual mean square","r-squared", "adjusted r-squared","AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant") 
tm=tmodel(m, data)
l1=list(l,con,tm,as.numeric(re1),re2[,1])
names(l1)=c("Schnute model","Parameters, standard error and confidence intervals","Whole model test","Residuals of schnute model", "Residuals standardized of schnute model") 
return(l1)
    }
    
    

fff2=function(m,m1,m2, name="Linear Model"){
name=name
# m model (random)
c1=fixef(m)[[1]]
c2=fixef(m)[[2]]
p=summary(m)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
A1=AIC(m)
B1=BIC(m)
r1=R2mm(m)
ar1=R3mm(m)
fr1=fr(m)
estimates1=c(c1,c2,p1,p2,r1,ar1,A1,B1,fr1); estimates1=round(estimates1,digits)
# m1 model (a random)
c1=fixef(m1)[[1]]
c2=fixef(m1)[[2]]
p=summary(m1)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
A1=AIC(m1)
B1=BIC(m1)
r1=R2mm(m1)
ar1=R3mm(m1)
fr1=fr(m1)
estimates2=c(c1,c2,p1,p2,r1,ar1,A1,B1,fr1); estimates2=round(estimates2,digits)
# m2 model (b random)
c1=fixef(m2)[[1]]
c2=fixef(m2)[[2]]
p=summary(m2)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
A1=AIC(m2)
B1=BIC(m2)
r1=R2mm(m2)
ar1=R3mm(m2)
fr1=fr(m2)
con=fcon(m)
con1=fcon(m1)
con2=fcon(m2)
fw=fwm(m);fw1=fwm(m1);fw2=fwm(m2)
estimates3=c(c1,c2,p1,p2,r1,ar1,A1,B1,fr1); estimates3=round(estimates3,digits)
resp1=data.frame(estimates1,estimates2,estimates3);names(resp1)=c("model random","a random","b random")
rownames(resp1)=c("coefficient a","coefficient b","p-value t.test for a","p-value t.test for b", "r-squared","adjusted r-squared", "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
lm1=as.data.frame(coef(m))
lm2=as.data.frame(coef(m1))
lm3=as.data.frame(coef(m2))
rr1=resid(m);rr2=resid(m1);rr3=resid(m2)
rr1p=scale(rr1)[,1];rr2p=scale(rr2)[,1];rr3p=scale(rr3)[,1]
da1=data.frame(rr1,rr2,rr3);da1=round(da1,digits);da2=data.frame(rr1p,rr2p,rr3p)
da2=round(da2,digits);names(da1)=c("model random","a random","b random")
names(da2)=c("model random","a random","b random")
rl=list(resp1,lm1,lm2,lm3,con,con1,con2,fw,fw1,fw2,da1,da2); names(rl)=c(name, "Coefficientes (model random)","Coefficientes (a random)","Coefficientes (b random)","Parameters, standard error and confidence intervals (model random)","Parameters, standard error and confidence intervals (a random)","Parameters, standard error and confidence intervals (b random)" ,"Anova (model random)","Anova (a random)","Anova (b random)","Residuals","Residuals (standardized)") 
return(rl)
}

fff3=function(m,m1,m2,m3, name="Quadratic Model"){
name=name
# m model (random)
c1=fixef(m)[[1]]
c2=fixef(m)[[2]]
c3=fixef(m)[[3]]
p=summary(m)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
A1=AIC(m)
B1=BIC(m)
r1=R2mm(m)
ar1=R3mm(m)
fr1=fr(m)
estimates1=c(c1,c2,c3,p1,p2,p3,r1,ar1,A1,B1,fr1); estimates1=round(estimates1,digits)
# m1 a (random)
c1=fixef(m1)[[1]]
c2=fixef(m1)[[2]]
c3=fixef(m1)[[3]]
p=summary(m1)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
A1=AIC(m1)
B1=BIC(m1)
r1=R2mm(m1)
ar1=R3mm(m1)
fr1=fr(m1)
estimates2=c(c1,c2,c3,p1,p2,p3,r1,ar1,A1,B1,fr1); estimates2=round(estimates2,digits)
# m2 b (random)
c1=fixef(m2)[[1]]
c2=fixef(m2)[[2]]
c3=fixef(m2)[[3]]
p=summary(m2)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
A1=AIC(m2)
B1=BIC(m2)
r1=R2mm(m2)
ar1=R3mm(m2)
fr1=fr(m2)
estimates3=c(c1,c2,c3,p1,p2,p3,r1,ar1,A1,B1,fr1); estimates3=round(estimates3,digits)
# m3 c (random)
c1=fixef(m3)[[1]]
c2=fixef(m3)[[2]]
c3=fixef(m3)[[3]]
p=summary(m3)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
A1=AIC(m3)
B1=BIC(m3)
r1=R2mm(m3)
ar1=R3mm(m3)
fr1=fr(m3)
estimates4=c(c1,c2,c3,p1,p2,p3,r1,ar1,A1,B1,fr1); estimates4=round(estimates4,digits)
resp1=data.frame(estimates1,estimates2,estimates3, estimates4);names(resp1)=c("model random","a random","b random", "c random")
rownames(resp1)=c("coefficient a","coefficient b","coefficient c","p-value t.test for a","p-value t.test for b","p-value t.test for c", "r-squared","adjusted r-squared", "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
lm1=as.data.frame(coef(m))
lm2=as.data.frame(coef(m1))
lm3=as.data.frame(coef(m2))
lm4=as.data.frame(coef(m3))
con=fcon(m);con1=fcon(m1);con2=fcon(m2); con3=fcon(m3)
fw=fwm(m);fw1=fwm(m1);fw2=fwm(m2); fw3=fwm(m3)
rr1=resid(m);rr2=resid(m1);rr3=resid(m2);rr4=resid(m3)
rr1p=scale(rr1)[,1];rr2p=scale(rr2)[,1];rr3p=scale(rr3)[,1];rr4p=scale(rr4)[,1]
da1=data.frame(rr1,rr2,rr3,rr4);da1=round(da1,digits);da2=data.frame(rr1p,rr2p,rr3p,rr4p)
da2=round(da2,digits);names(da1)=c("model random","a random","b random","c random")
names(da2)=c("model random","a random","b random","c random")
rl=list(resp1,lm1,lm2,lm3,lm4,con,con1,con2,con3,fw,fw1,fw2,fw3,da1,da2); names(rl)=c(name, "Coefficientes (model random)","Coefficientes (a random)","Coefficientes (b random)","Coefficientes (c random)","Parameters, standard error and confidence intervals (model random)","Parameters, standard error and confidence intervals (a random)" , "Parameters, standard error and confidence intervals (b random)","Parameters, standard error and confidence intervals (c random)","Anova (model random)","Anova (a random)","Anova (b random)","Anova (c random)","Residuals","Residuals (standardized)") 
return(rl)
}

fff4=function(m,m1,m2,m3,m4, name="Bi Linear Model"){
name=name
# m model (random)
c1=fixef(m)[[1]]
c2=fixef(m)[[2]]
c3=fixef(m)[[3]]
c4=fixef(m)[[4]]
p=summary(m)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
A1=AIC(m)
B1=BIC(m)
r1=R2mm(m)
ar1=R3mm(m)
fr1=fr(m)
estimates1=c(c1,c2,c3,c4,p1,p2,p3,p4,r1,ar1,A1,B1,fr1); estimates1=round(estimates1,digits)
# m1 a (random)
c1=fixef(m1)[[1]]
c2=fixef(m1)[[2]]
c3=fixef(m1)[[3]]
c4=fixef(m1)[[4]]
p=summary(m1)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
A1=AIC(m1)
B1=BIC(m1)
r1=R2mm(m1)
ar1=R3mm(m1)
fr1=fr(m1)
estimates2=c(c1,c2,c3,c4,p1,p2,p3,p4,r1,ar1,A1,B1,fr1); estimates2=round(estimates2,digits)
# m2 b (random)
c1=fixef(m2)[[1]]
c2=fixef(m2)[[2]]
c3=fixef(m2)[[3]]
c4=fixef(m2)[[4]]
p=summary(m2)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
A1=AIC(m2)
B1=BIC(m2)
r1=R2mm(m2)
ar1=R3mm(m2)
fr1=fr(m2)
estimates3=c(c1,c2,c3,c4,p1,p2,p3,p4,r1,ar1,A1,B1,fr1); estimates3=round(estimates3,digits)
# m3 c (random)
c1=fixef(m3)[[1]]
c2=fixef(m3)[[2]]
c3=fixef(m3)[[3]]
c4=fixef(m3)[[4]]
p=summary(m3)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
A1=AIC(m3)
B1=BIC(m3)
r1=R2mm(m3)
ar1=R3mm(m3)
fr1=fr(m3)
estimates4=c(c1,c2,c3,c4,p1,p2,p3,p4,r1,ar1,A1,B1,fr1); estimates4=round(estimates4,digits)
# m4 d (random)
c1=fixef(m4)[[1]]
c2=fixef(m4)[[2]]
c3=fixef(m4)[[3]]
c4=fixef(m4)[[4]]
p=summary(m4)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
A1=AIC(m4)
B1=BIC(m4)
r1=R2mm(m4)
ar1=R3mm(m4)
fr1=fr(m4)
estimates5=c(c1,c2,c3,c4,p1,p2,p3,p4,r1,ar1,A1,B1,fr1); estimates5=round(estimates5,digits)
resp1=data.frame(estimates1,estimates2,estimates3, estimates4, estimates5);names(resp1)=c("model random","a random","b random", "c random", "d random")
rownames(resp1)=c("coefficient a","coefficient b","coefficient c","coefficient d","p-value t.test for a","p-value t.test for b","p-value t.test for c","p-value t.test for d", "r-squared","adjusted r-squared", "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
lm1=as.data.frame(coef(m))
lm2=as.data.frame(coef(m1))
lm3=as.data.frame(coef(m2))
lm4=as.data.frame(coef(m3))
lm5=as.data.frame(coef(m4))
con=fcon(m);con1=fcon(m1);con2=fcon(m2); con3=fcon(m3);con4=fcon(m4)
fw=fwm(m);fw1=fwm(m1);fw2=fwm(m2); fw3=fwm(m3);fw4=fwm(m4)
rr1=resid(m);rr2=resid(m1);rr3=resid(m2);rr4=resid(m3);rr5=resid(m4)
rr1p=scale(rr1)[,1];rr2p=scale(rr2)[,1];rr3p=scale(rr3)[,1];rr4p=scale(rr4)[,1];rr5p=scale(rr5)[,1]
da1=data.frame(rr1,rr2,rr3,rr4,rr5);da1=round(da1,digits);da2=data.frame(rr1p,rr2p,rr3p,rr4p,rr5p)
da2=round(da2,digits);names(da1)=c("model random","a random","b random","c random","d random")
names(da2)=c("model random","a random","b random","c random","d random")
rl=list(resp1,lm1,lm2,lm3,lm4,lm5,con,con1,con2,con3,con4,fw,fw1,fw2,fw3,fw4,da1,da2); names(rl)=c(name, "Coefficientes (model random)","Coefficientes (a random)","Coefficientes (b random)","Coefficientes (c random)","Coefficientes (d random)", "Parameters, standard error and confidence intervals (model random)","Parameters, standard error and confidence intervals (a random)" , "Parameters, standard error and confidence intervals (b random)","Parameters, standard error and confidence intervals (c random)","Parameters, standard error and confidence intervals (d random)","Anova (model random)","Anova (a random)","Anova (b random)","Anova (c random)","Anova (d random)","Residuals","Residuals (standardized)") 
return(rl)
}

fff5=function(m,m1,m2,m3,m4,m5, name="Logistic Bi Compartmental Model"){
name=name
# m model (random)
c1=fixef(m)[[1]]
c2=fixef(m)[[2]]
c3=fixef(m)[[3]]
c4=fixef(m)[[4]]
c5=fixef(m)[[5]]
p=summary(m)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m)
B1=BIC(m)
r1=R2mm(m)
ar1=R3mm(m)
fr1=fr(m)
estimates1=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates1=round(estimates1,digits)
# m1 a (random)
c1=fixef(m1)[[1]]
c2=fixef(m1)[[2]]
c3=fixef(m1)[[3]]
c4=fixef(m1)[[4]]
c5=fixef(m1)[[5]]
p=summary(m1)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m1)
B1=BIC(m1)
r1=R2mm(m1)
ar1=R3mm(m1)
fr1=fr(m1)
estimates2=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates2=round(estimates2,digits)
# m2 b (random)
c1=fixef(m2)[[1]]
c2=fixef(m2)[[2]]
c3=fixef(m2)[[3]]
c4=fixef(m2)[[4]]
c5=fixef(m2)[[5]]
p=summary(m2)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m2)
B1=BIC(m2)
r1=R2mm(m2)
ar1=R3mm(m2)
fr1=fr(m2)
estimates3=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates3=round(estimates3,digits)
# m3 c (random)
c1=fixef(m3)[[1]]
c2=fixef(m3)[[2]]
c3=fixef(m3)[[3]]
c4=fixef(m3)[[4]]
c5=fixef(m3)[[5]]
p=summary(m3)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m3)
B1=BIC(m3)
r1=R2mm(m3)
ar1=R3mm(m3)
fr1=fr(m3)
estimates4=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates4=round(estimates4,digits)
# m4 d (random)
c1=fixef(m4)[[1]]
c2=fixef(m4)[[2]]
c3=fixef(m4)[[3]]
c4=fixef(m4)[[4]]
c5=fixef(m4)[[5]]
p=summary(m4)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m4)
B1=BIC(m4)
r1=R2mm(m4)
ar1=R3mm(m4)
fr1=fr(m4)
estimates5=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates5=round(estimates5,digits)
# m5 e (random)
c1=fixef(m5)[[1]]
c2=fixef(m5)[[2]]
c3=fixef(m5)[[3]]
c4=fixef(m5)[[4]]
c5=fixef(m5)[[5]]
p=summary(m5)$tTable[,5]
p1=p[[1]]
p2=p[[2]]
p3=p[[3]]
p4=p[[4]]
p5=p[[5]]
A1=AIC(m5)
B1=BIC(m5)
r1=R2mm(m5)
ar1=R3mm(m5)
fr1=fr(m5)
estimates6=c(c1,c2,c3,c4,c5,p1,p2,p3,p4,p5,r1,ar1,A1,B1,fr1); estimates6=round(estimates6,digits)
resp1=data.frame(estimates1,estimates2,estimates3, estimates4, estimates5, estimates6);names(resp1)=c("model random","a random","b random", "c random", "d random", "e random")
rownames(resp1)=c("coefficient a","coefficient b","coefficient c","coefficient d","coefficient e","p-value t.test for a","p-value t.test for b","p-value t.test for c","p-value t.test for d","p-value t.test for e", "r-squared","adjusted r-squared", "AIC", "BIC", "p.value Shapiro-Wilk test","coefficient of variation (%)", "first value most discrepant","second value most discrepant","third value most discrepant")
lm1=as.data.frame(coef(m))
lm2=as.data.frame(coef(m1))
lm3=as.data.frame(coef(m2))
lm4=as.data.frame(coef(m3))
lm5=as.data.frame(coef(m4))
lm6=as.data.frame(coef(m5))
con=fcon(m);con1=fcon(m1);con2=fcon(m2); con3=fcon(m3);con4=fcon(m4); con5=fcon(m5)
fw=fwm(m);fw1=fwm(m1);fw2=fwm(m2); fw3=fwm(m3);fw4=fwm(m4);fw5=fwm(m5)
rr1=resid(m);rr2=resid(m1);rr3=resid(m2);rr4=resid(m3);rr5=resid(m4);rr6=resid(m5)
rr1p=scale(rr1)[,1];rr2p=scale(rr2)[,1];rr3p=scale(rr3)[,1];rr4p=scale(rr4)[,1]
rr5p=scale(rr5)[,1];rr6p=scale(rr6)[,1]
da1=data.frame(rr1,rr2,rr3,rr4,rr5,rr6);da1=round(da1,digits);da2=data.frame(rr1p,rr2p,rr3p,rr4p,rr5p,rr6p)
da2=round(da2,digits);names(da1)=c("model random","a random","b random","c random","d random","e random")
names(da2)=c("model random","a random","b random","c random","d random","e random")
rl=list(resp1,lm1,lm2,lm3,lm4,lm5,lm6,con,con1,con2,con3,con4,con5,fw,fw1,fw2,fw3,fw4,fw5,da1,da2); names(rl)=c(name, "Coefficientes (model random)","Coefficientes (a random)","Coefficientes (b random)","Coefficientes (c random)","Coefficientes (d random)","Coefficientes (e random)", "Parameters, standard error and confidence intervals (model random)","Parameters, standard error and confidence intervals (a random)" , "Parameters, standard error and confidence intervals (b random)","Parameters, standard error and confidence intervals (c random)","Parameters, standard error and confidence intervals (d random)","Parameters, standard error and confidence intervals (e random)","Anova (model random)","Anova (a random)","Anova (b random)","Anova (c random)","Anova (d random)","Anova (e random)","Residuals","Residuals (standardized)") 
return(rl)
}

# functions for mixed models
lf=list(
ff1=function(x, a,b){a+b*x}, # linear
ff2=function(x, a,b,c){a+b*x+c*x^2}, # quadratic
ff3=function(x, a,b,c){a+b*(x-c)*(x<=c)}, # plateau l
ff4=function(x, a,b,c){(a+b*x+c*I(x^2))*(x<=-0.5*b/c)+(a + I(-b^2/(4*c)))*(x>-0.5*b/c)}, 
# plateau q
ff5=function(x,a,b,c,d){ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)},# bi linear
ff6=function(x, a,b){a*exp(b*x)}, # exponential
ff7=function(x, a,b,c){a*(1+b*(exp(-c*x)))^-1}, # logistic
ff8=function(x, a,b,c){a*(1-b*(exp(-c*x)))^3}, # Von
ff9=function(x, a,b,c){a*(1-b*(exp(-c*x)))}, # brody
ff10=function(x, a,b,c){a*exp(-b*exp(-c*x))}, # gompertz
ff11=function(x, a,b,c){(a*x^b)*exp(-c*x)}, # lactation
ff12=function(x, a,b,c){a + b * (1 - exp(-c * x))}, # orskov
ff13=function(x, a,b,c,d,e){(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e))))}, # bi logistic
ff14=function(x, a,b){a*(b*x)}, # exponential (allometric)
ff15=function(x, a,b,c,d){a+b*x+c*x^2+d*x^3}, # cubic
ff16=function(x, a,b,c){a/(1+b*(exp(-c*x)))^d} # richards
)

f1m=function(data){
names(data)=c("x","blocks","y")
data$blocks=as.factor(data$blocks)

m=lme(y ~ x, random = ~ x| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m1=lme(y ~ x, random = ~ 1| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m2=lme(y ~ x, random = ~ -1+x| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 
fff2(m,m1,m2,name="Linear Model")
}

# quadratic
f2m=function(data){
names(data)=c("x","blocks","y")
data$blocks=as.factor(data$blocks)

m=lme(y ~ x+I(x^2), random = ~ x+I(x^2)| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m1=lme(y ~ x+I(x^2), random = ~ 1| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m2=lme(y ~ x+I(x^2), random = ~ -1+x| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m3=lme(y ~ x+I(x^2), random = ~ -1+I(x^2)| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 
gr=fff3(m,m1,m2,m3,name="Quadratic Model")
return(gr)
}

# linear plateau
f3m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
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
l1=nlsList(y~a+b*(x-c)*(x<=c)|blocks, data=data, start=list(a=a11,b=b11,c=c11),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m1=nlme(l1, fixed=a+b+c~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m2=nlme(l1, fixed=a+b+c~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m3=nlme(l1, fixed=a+b+c~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Linear Plateau Model")
return(gr)
}


# quadratic plateau
f4m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
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
l1=nlsList(y~(a+b*x+c*I(x^2))*(x<=-0.5*b/c)+(a + I(-b^2/(4*c)))*(x>-0.5*b/c)|blocks, data=data, start=list(a=a11,b=b11,c=c11),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m1=nlme(l1, fixed=a+b+c~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m2=nlme(l1, fixed=a+b+c~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m3=nlme(l1, fixed=a+b+c~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Quadratic Plateau Model")
return(gr)
}

# bi linear 
f5m=function(data){
data=data
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3],d=s[4]),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c+d~x, random=a+b+c+d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m1=nlme(l1, fixed=a+b+c+d~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m2=nlme(l1, fixed=a+b+c+d~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m3=nlme(l1, fixed=a+b+c+d~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m4=nlme(l1, fixed=a+b+c+d~x, random=d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff4(m,m1,m2,m3,m4,name="Broken Line Model")
return(gr)
}

# exponential
f6m=function(data){
names(data)=c("x","blocks","y")
data$blocks=as.factor(data$blocks)
l = nlsList(y~a*exp(b*x)|blocks, data=data,
start=list(a=s[1], b=s[2]),control = nls.control(maxiter = 6000),na.action = na.omit)

m <- nlme(l, fixed=a+b~x, random=a+b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m1 <- nlme(l, fixed=a+b~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m2 <- nlme(l, fixed=a+b~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff2(m,m1,m2,name="Exponential Model")
return(gr)
}

# logistic
f7m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~a*(1+b*(exp(-c*x)))^-1|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m1=nlme(l1, fixed=a+b+c~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m2=nlme(l1, fixed=a+b+c~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m3=nlme(l1, fixed=a+b+c~x, random=c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Logistic Model")
return(gr)
}

# Van B
f8m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~a*(1-b*(exp(-c*x)))^3|blocks, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c~x, random=c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Van Bertalanffy Model")
return(gr)
}

# brody
f9m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~a*(1-b*(exp(-c*x)))|blocks, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c~x, random=c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Brody Model")
return(gr)
}

# gompertz
f10m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~a*exp(-b*exp(-c*x))|blocks, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c~x, random=c~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Gompertz Model")
return(gr)
}

# lactation
f11m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~(a*x^b)*exp(-c*x)|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3]),control = nls.control(maxiter = 6000),na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c~x, random=a~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c~x, random=b~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c~x, random=c~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Lactation Model")
return(gr)
}

# orskov
f12m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~a + b * (1 - exp(-c * x))|blocks, data=data, start=list(a=20,b=60,c=0.05),control = nls.control(maxiter = 6000),na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c~x, random=a+b+c~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c~x, random=a~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c~x, random=b~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c~x, random=c~1, data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff3(m,m1,m2,m3,name="Ruminal Degradation Model")
return(gr)
}

# logistic bi-compartmental
f13m=function(data){
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e))))|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3],d=s[4],e=s[5]),control = nls.control(maxiter = 6000),na.action = na.omit)
# 
m=nlme(l1, fixed=a+b+c+d+e~x, random=a+b+c+d+e~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# a 
m1=nlme(l1, fixed=a+b+c+d+e~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# b 
m2=nlme(l1, fixed=a+b+c+d+e~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# c 
m3=nlme(l1, fixed=a+b+c+d+e~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# d 
m4=nlme(l1, fixed=a+b+c+d+e~x, random=d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
# e 
m5=nlme(l1, fixed=a+b+c+d+e~x, random=e~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff5(m,m1,m2,m3,m4,m5,name="Logistic Bi-Compartmental Model")
return(gr)
}

# exponential (allometric)
f14m=function(data){
names(data)=c("x","blocks","y")
data$blocks=as.factor(data$blocks)
l = nlsList(y~a*(x^b)|blocks, data=data,
start=list(a=s[1], b=s[2]),control = nls.control(maxiter = 6000),na.action = na.omit)

m <- nlme(l, fixed=a+b~x, random=a+b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m1 <- nlme(l, fixed=a+b~x, random=a~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)

m2 <- nlme(l, fixed=a+b~x, random=b~1,control = nls.control(maxiter = 6000),data=data,na.action = na.omit)
gr=fff2(m,m1,m2,name="Exponential Model (allometric model)")
return(gr)
}

# cubic
f15m=function(data){
names(data)=c("x","blocks","y")
data$blocks=as.factor(data$blocks)

m=lme(y ~ x+I(x^2)+I(x^3), random = ~ x+I(x^2)+I(x^3)| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m1=lme(y ~ x+I(x^2)+I(x^3), random = ~ 1| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m2=lme(y ~ x+I(x^2)+I(x^3), random = ~ -1+x| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m3=lme(y ~ x+I(x^2)+I(x^3), random = ~ -1+I(x^2)| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 

m4=lme(y ~ x+I(x^2)+I(x^3), random = ~ -1+I(x^3)| blocks, data=data,method="REML",control = lmeControl(maxIter = 6000, msMaxIter = 6000, niterEM = 2000, opt = "optim"),na.action = na.omit) 


gr=fff4(m,m1,m2,m3,m4, name="Cubic Model")
return(gr)
}

# richards
f16m=function(data){
data=data
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~(a/(1+b*(exp(-c*x)))^d)|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3],d=s[4]),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c+d~x, random=a+b+c+d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m1=nlme(l1, fixed=a+b+c+d~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m2=nlme(l1, fixed=a+b+c+d~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m3=nlme(l1, fixed=a+b+c+d~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m4=nlme(l1, fixed=a+b+c+d~x, random=d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff4(m,m1,m2,m3,m4,name="Richards Model")
return(gr)
}

# schnute
f17m=function(data){
data=data
names(data)=c("x","blocks","y")
data[,2]=factor(data[,2])
l1=nlsList(y~(a^d+ ((b^d)-(a^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d)|blocks, data=data, start=list(a=s[1],b=s[2],c=s[3],d=s[4]),control = nls.control(maxiter = 6000),na.action = na.omit)

m=nlme(l1, fixed=a+b+c+d~x, random=a+b+c+d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m1=nlme(l1, fixed=a+b+c+d~x, random=a~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m2=nlme(l1, fixed=a+b+c+d~x, random=b~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m3=nlme(l1, fixed=a+b+c+d~x, random=c~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)

m4=nlme(l1, fixed=a+b+c+d~x, random=d~1,data=data,control = nls.control(maxiter = 6000),na.action = na.omit)
gr=fff4(m,m1,m2,m3,m4,name="Schnute Model")
return(gr)
}
   
    reg=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14, f15, f16, f17)
    funreg=list(f1m,f2m,f3m,f4m,f5m,f6m,f7m,f8m,f9m,f10m,f11m,f12m,f13m,f14m, f15m, f16m, f17m)    
i=ifelse(mixed==FALSE,1,2)
lf=list(reg,funreg)
gg=lf[[i]]
FUN=gg[[model]]

pf=list(opt1,opt2)
pf=pf[[i]]
pf=pf(data)
rep=lapply(pf, FUN)
ln=list(1,c(1,2))
ln=ln[[i]]
names(rep)= names(data)[-ln]
return(rep)   
}
