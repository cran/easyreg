bl<-function (data, model=1, alpha=0.05, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20, xlim="default.x",ylim="default.y",...) 
{
data=na.omit(data);model=model; xlab = xlab; ylab=ylab;position = position; digits = digits;mean = mean
sd=sd; legend = legend; lty=lty; col=col; pch=pch

tm2=function(npar,xlin,ixx,res){
by=(xlin%*%data$y)
par=ixx%*%by
sqt=(t(par)%*%by)-((sum(data$y)^2)/length(data$y))
sqr=(t(par[c(1:npar),])%*%by[c(1:npar),])-((sum(data$y)^2)/length(data$y))
sqb=sqt-sqr
sqe=sum(res^2)
gle=(length(res)-(npar+(nlevels(data$blocks)-1)))
sqt=var(data$y)*(length(data$y)-1)
sq=c(sqr, sqb=sqt-sqe-sqr,sqe, sqt)
df=c(npar-1,nlevels(data$blocks)-1,gle, length(data$y)-1)
qm=sq/df
fc=qm[1]/qm[3]
p=pf(fc,2,gle, lower.tail=FALSE)
f=c(round(fc,2),"-","-","-")
p=if(p<0.001) "<0.001" else p=round(p,4)
p=c(p,"-","-","-")
tab=data.frame(df,round(sq,4),round(qm,4),f,p)
source=c("regression", "blocks","residuals", "total")
colnames(tab)=c("df","sum of squares", "mean squares", "Fcal", "p-value")
rownames(tab)=source
return(tab)
}


fcon<-function(m){
s=diag(vcov(m))^0.5
i=if(class(m)=="nls") 1 else 2
gl=list(summary(m)$df[2],m$df.residual)
gl=gl[[i]]
t=abs(qt((alpha/2),gl))
c=as.numeric(coef(m))
di=t*s
c1=c-di
c2=c+di
l=letters[1:length(coef(m))]
r=data.frame(l,c,s,c1,c2)	
p1=paste("inf",(alpha/2)*100,sep=".")
p2=paste("sup",(1-(alpha/2))*100,sep=".")
names(r)=c("parameter", "estimate", "standard error", p1,p2) 
r=data.frame(Parameters=r[,1],round(r[,-1], digits))
rownames(r)=NULL
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


fcon2<-function(gl,cof,s){
t=abs(qt((alpha/2),gl))
c=as.numeric(cof)
di=t*s
c1=c-di
c2=c+di
l=letters[1:length(c)]
r=data.frame(l,c,s,c1,c2)	
p1=paste("inf",(alpha/2)*100,sep=".")
p2=paste("sup",(1-(alpha/2))*100,sep=".")
names(r)=c("parameter", "estimate", "standard error", p1,p2) 
r=data.frame(Parameters=r[,1],round(r[,-1], digits))
rownames(r)=NULL
return(r)
}

    	minx = min(data[, 1], na.rm=TRUE) - sd(data[, 1],na.rm=TRUE)/2
        maxx = max(data[, 1],na.rm=TRUE) + sd(data[, 1],na.rm=TRUE)/2
        miny = min(if(ncol(data)==2)data[, 2]else data[,3],na.rm=TRUE) - sd(if(ncol(data)==2)data[, 2]else data[,3],na.rm=TRUE)/2
        maxy = max(if(ncol(data)==2)data[, 2]else data[,3],na.rm=TRUE) + sd(if(ncol(data)==2)data[, 2]else data[,3],na.rm=TRUE)/2
        
xli=ifelse(xlim=="default.x",1,2)
xli=xli[1]
xll=list(c(minx,maxx),xlim)
xl=xll[[xli]]

yli=ifelse(ylim=="default.y",1,2)
yli=yli[1]
yll=list(c(miny,maxy),ylim)
yl=yll[[yli]]


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

tmf<-function(res, data, df)
{
sqt=var(data$y)*(length(data$y)-1)
sqres=sum(res^2)
sqr=sqt-sqres
df1=df
df2=length(data$y)-1
df3=df2-df1
df=c(df1,df3,df2)
sq=c(sqr,sqres,sqt)
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


tmodel<-function(m, data){
m0=lm(y~1, data=data)
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


fol=function(data){
yy=ifelse(ncol(data)==2, 2,3)
# ornenar por x
data=data[order(data[,1]),]
# criar uma coluna fator de x
data=data.frame(fac=as.factor(data[,1]),data)
# slice de x
l=1:nlevels(data[,1])
fs1=function(i){s=subset(data, data[,1]==data[,1][i]); return(s)}
s=lapply(l, fs1)
# retirar a coluna do fator
ifac=1:length(s)
ffac=function(i){return(s[[i]][,-1])}
lfac=lapply(ifac,ffac)
# juntar uma a uma
i=1:length(lfac)
n=length(i)
i=i[-c(1,n,n-1)]
f1=function(i){sl=lfac[c(1:i)]; return(sl)}
dl=lapply(i, f1)
j=1:length(dl)
fc1=function(i){u=unlist(dl[[i]]);u=length(u)/yy;return(u)}
dll=lapply(j, fc1)
fl1=function(i){u=data[c(1:i),-1];return(u)}
dl=lapply(dll,fl1)
i=3:(n-1)
f2=function(i){sl=lfac[c(i:n)]; return(sl)}
dp=lapply(i, f2)
fc2=function(i){u=unlist(dp[[i]]);u=length(u)/yy;return(u)}
dlp=lapply(j, fc2)
fl1=function(i){u=data[-c(1:i),-1];return(u)}
dp=lapply(dll,fl1)
data=data[,-1]
l=list(dl,dp, data); return(l)
}


fop=function(data)
{
yy=ifelse(ncol(data)==2, 2,3)
# ornenar por x
data=data[order(data[,1]),]
# criar uma coluna fator de x
data=data.frame(fac=as.factor(data[,1]),data)
# slice de x
l=1:nlevels(data$fac)
fs=function(i){s=subset(data, data[,1]==data$fac[i]); return(s)}
s=lapply(l, fs)
# retirar a coluna do fator
ifac=1:length(s)
ffac=function(i){return(s[[i]][,-1])}
lfac=lapply(ifac,ffac)
# juntar uma a uma
i=1:length(lfac)
n=length(i)
i=i[-c(1,2,n,n-1)]
f1=function(i){sl=lfac[c(1:i)]; return(sl)}
dl=lapply(i, f1)
j=1:length(dl)
fc1=function(i){u=unlist(dl[[i]]);u=length(u)/yy;return(u)}
dll=lapply(j, fc1)
fl1=function(i){u=data[c(1:i),-1];return(u)}
dl=lapply(dll,fl1)
i=4:(n-1)
f2=function(i){sl=lfac[c(i:n)]; return(sl)}
dp=lapply(i, f2)
fc2=function(i){u=unlist(dp[[i]]);u=length(u)/yy;return(u)}
dlp=lapply(j, fc2)
fl1=function(i){u=data[-c(1:i),-1];return(u)}
dp=lapply(dll,fl1)
data=data[,-1]
l=list(dl,dp, data); return(l)
}



means = function(data) {
	data=if(ncol(data)==3) data[,-2] else data
        t = as.factor(data[, 1])
        d = data.frame(t, data[, -1])
        s = split(data.frame(d[, -1]), d$t)
        r = lapply(s, colMeans, na.rm = TRUE)
        r = lapply(r, round, 2)
	fg=function(i){ss=s[[i]][,1]; return(ss)}
	i=1:nlevels(t)
	r2=lapply(i, fg)
	r2 = lapply(r2, sd, na.rm = TRUE)
	r2=unlist(r2)
        rr = t(data.frame(r))
        rr = data.frame(rr)
        rownames(rr) = NULL
        treat = levels(t)
        rr = data.frame(treat, rr, r2)
        colnames(rr) = c("x","y","sd" )
        return(rr)
    }

bl1=function (data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...)
{
	sdd=sd   
	names(data)=c("x","y")
	mini=mean(data$x, na.rm=TRUE)
	mini=(mini/100)
	knot=seq(from=min(data$x), to=max(data$x), by=mini)
	fm=function(knot){m1=lm(y~x+(pmax(0,x-knot, na.rm=TRUE)), data=data);r=sum(resid(m1)^2);return(r)}
	s1=sapply(knot, fm)
	s=data.frame(knot,s1)
	s=s[order(s$s1),]
	knot=s[1,1]
	m=lm(y~x+(pmax(0,x-knot, na.rm=TRUE)), data=data)
	a=coef(m)[1]
	b=coef(m)[2]
	a=round(a,digits)
	b=round(b,digits)
	f1=function(x){a+b*x}
	d=(coef(m)[2]+coef(m)[3])
	resf1=f1(knot)
	c=resf1-(d*knot)
	c=round(c,digits)
	d=round(d,digits)
	f2=function(x){c+d*x}
	cofs=summary(m)[[4]]
	cofs=cofs[,c(1,4)]
	colnames(cofs)=c("coefficientes","p_values")
	rownames(cofs)=NULL;cofs=as.data.frame(cofs)
	l1=c(a,b);l2=c(c,d)
	l=data.frame(l1,l2)
	names(l)=c("first line","second line")
	rownames(l)=c("a","b")
	l=t(l)
	k=knot
	r1=summary(m)[[8]]
	r2=summary(m)[[9]]
	res=resid(m)
	shap=shapiro.test(resid(m))
	shap=shap$p.value
	sres=scale(res)
	sres=sres[,1]
	tm=tmodel(m,data)
	fc=fcon(m)
	resp=list(round(cofs,digits), round(l,digits), knot, tm, fc,round(r1,digits),round(r2,digits),round(AIC(m),digits),round(BIC(m),digits),res, sres, round(shap,digits) )
	names(resp)=c("Coefficients", "Lines","Knot (Break-Point)","Whole model test","Parameters, standard error and confidence intervals", "R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")

		sin1 = ifelse(a > 0, "+", "")
		sin2 = ifelse(b > 0, "+", "")
		sin3 = ifelse(c > 0, "+", "")
		sin4 = ifelse(d > 0, "+", "")

	r1=round(r1,2)
	knot=round(knot,digits)

e5 = substitute(atop(y==a*sin2*b*x,y==c*sin4*d*x)*"  "*atop(x<knot,x>=knot)*"    "* R^2 * " = " * r1, list(a = a, b = b, c = c,d = d, r1 = r1,knot=knot, sin2 = sin2, sin4 = sin4))

t = list("top", "bottomright", "bottom", "bottomleft", 
             "left", "topleft", "topright", "right", "center")
    p = t[[position]]


mei=mean
 
datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=data[,2], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
    
  
#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)


fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
    }
    
ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}



bl2=function (data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...) 
{
 sdd=sd 
names(data) = c("x", "y")
fo=fol(data)
dl=fo[[1]]
dp=fo[[2]]
data=fo[[3]]

fml=function(data){m=lm(y~x, data=data); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
fmp=function(data){m=lm(y~1, data=data); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
om=1:length(m1)
mr=data.frame(unlist(m1), unlist(m2));names(mr)=c("l","p") 
s=rowSums(mr)
mr=data.frame(mr,s, om)
mr = mr[order(mr$s), ]
bm=mr[1,4]
bm2=mr[2,4]
fml=function(data){m=lm(y~x, data=data); return(m)}
fmp=function(data){m=lm(y~1, data=data); return(m)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m2)[[1]]
knot=(c-a)/b
bm=ifelse(knot>min(d2[,1]),bm2,bm)
fml=function(data){m=lm(y~x, data=data); return(m)}
fmp=function(data){m=lm(y~1, data=data); return(m)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m2)[[1]]
knot=(c-a)/b
      cofs = c(a,b,c); names(cofs)=c("a","b","plateau")
    res = c(as.numeric(resid(m1)),as.numeric(resid(m2)))
    shap = shapiro.test(res)
    shap = shap$p.value
    sres = scale(res)
    sres = sres[, 1]
    vd=sum(res^2)/(length(res)-3)
rp1=rep(1,length(d1$x))
rp2=rep(0,length(d2$x))
rp3=d1$x
rp4=rep(0,length(d2$x))
rp5=rep(0,length(d1$x))
rp6=rep(1,length(d2$x))
rr=c(rp1,rp2,rp3,rp4,rp5,rp6)
x=matrix(rr, ncol=3)
xlin=t(x)
xx=xlin%*%x
ixx=solve(xx)
vco=ixx*vd; sp=diag(vco)^0.5
ta=a/sp[1];tb=b/sp[2]; tc=c/sp[3]
pa=pt(abs(ta),(length(res)-3), lower.tail = FALSE)
pa=pa*2
pb=pt(abs(tb),(length(res)-3), lower.tail = FALSE)
pb=pb*2
pc=pt(abs(tc),(length(res)-3), lower.tail = FALSE)
pc=pc*2
lll=c(pa,pb,pc)
names(lll)=c("a","b","plateau")
sqt = sum((data$y-mean(data$y))^2)
sqr=sqt-sum(res^2)
r1 = sqr/sqt
gl=length(res)-1
p1 =(gl/((gl + 1)-(3))*(1-r1))
r2 =1-p1
pred=na.exclude(data$y)-res
vdc=sum(res^2)/(length(res))    
sd=vdc^0.5 
loglik=sum(log(dnorm(x=data$y, mean=pred, sd=sd)))
k=3 
aic=-2*loglik+(2*(k+1))
bic=-2*loglik+(log(length(res))*(k+1))
x.plateau = knot
y.plateau = c
  tm=tmf(res,data,2)
fc=fcon2(length(res)-3,cofs,sp)

    resp = list(round(cofs, digits), round(lll, digits), knot,c,tm, 
        fc,round(r1, digits), round(r2, digits), round(aic, digits), 
        round(bic, digits), res, sres, round(shap, digits))
    names(resp) = c("Coefficients", "p values for coefficients", "x plateau","y plateau", 
        "Whole model test","Parameters, standard error and confidence intervals", "R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", 
        "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")
    sin1 = ifelse(a > 0, "+", "")
    sin2 = ifelse(b > 0, "+", "")
    sin3 = ifelse(c > 0, "+", "")
    r1 = round(r1, 2)
    knot = round(knot, digits)
    e5 = substitute(atop(y == a * sin2 * b * x, y == c) * "  " * atop(x < knot, x >= knot) * "    " * 
        R^2 * " = " * r1, list(a = a, b = b, c = c, r1 = r1, 
        knot = knot, sin2 = sin2))
    t = list("top", "bottomright", "bottom", "bottomleft", "left", 
        "topleft", "topright", "right", "center")
    p = t[[position]]
    mei = mean


###########################################
    
    datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=data[,2], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
    

#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)
##############################################

f1l=function(x)a+b*x
f2p=function(x)c+0*x
    fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,xlab = xlab,ylab=ylab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
    }
    
ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}


bl3=function(data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...)
{
sdd=sd       
names(data) = c("x", "y")
fo=fop(data)
dl=fo[[1]]
dp=fo[[2]]
data=fo[[3]]
   
fmq=function(data){m=lm(y~x+I(x^2), data=data);r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
fmp=function(data){m=lm(y~1, data=data); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
m1=lapply(dl, fmq)
m2=lapply(dp, fmp)
om=1:length(m1)
mr=data.frame(unlist(m1), unlist(m2));names(mr)=c("l","p") 
s=rowSums(mr)
mr=data.frame(mr,s, om)
mr = mr[order(mr$s), ]
bm=mr[1,4]
mq=function(data){m=lm(y~x+I(x^2), data=data); return(m)}
m1=lapply(dl, mq)
mp=function(data){m=lm(y~1, data=data); return(m)}
m2=lapply(dp, mp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m1)[[3]]
k=coef(m2)[[1]]
delta=sqrt((b^2)-(4*c*a)+(4*c*k))
knot1=((-1*b)+delta)/(2*c)
knot2=((-1*b)-delta)/(2*c)
k1=abs(knot1)
k2=abs(knot2)
kk=c(k1,k2); knot=sort(kk)[1]
    cofs = c(a,b,c,k); names(cofs)=c("a","b","c","plateau")
    res = c(as.numeric(resid(m1)[1:length(d1$x)]),as.numeric(resid(m2)))
    shap = shapiro.test(res)
    shap = shap$p.value
    sres = scale(res)
    sres = sres[, 1]
    vd=sum(res^2)/(length(res)-4)
rp1=rep(1,length(d1$x))
rp2=rep(0,length(d2$x))
rp3=d1$x
rp4=rep(0,length(d2$x))
rp5=d1$x^2
rp6=rep(0,length(d2$x))
rp7=rep(0,length(d1$x))
rp8=rep(1,length(d2$x))
rr=c(rp1,rp2,rp3,rp4,rp5,rp6,rp7,rp8)
x=matrix(rr, ncol=4)
xlin=t(x)
xx=xlin%*%x
ixx=solve(xx)
vco=ixx*vd; sp=diag(vco)^0.5
ta=a/sp[1];tb=b/sp[2]; tc=c/sp[3]; tk=k/sp[4]
pa=pt(abs(ta),(length(res)-4), lower.tail = FALSE)
pa=pa*2
pb=pt(abs(tb),(length(res)-4), lower.tail = FALSE)
pb=pb*2
pc=pt(abs(tc),(length(res)-4), lower.tail = FALSE)
pc=pc*2
pk=pt(abs(tk),(length(res)-4), lower.tail = FALSE)
pk=pk*2
lll=c(pa,pb,pc, pk)
names(lll)=c("a","b","c","plateau")
sqt = sum((data$y-mean(data$y))^2)
sqr=sqt-sum(res^2)
r1 = sqr/sqt
gl=length(res)-1
p1 =(gl/((gl + 1)-(4))*(1-r1))
r2 =1-p1
pred=data$y-res
vdc=sum(res^2)/(length(res))    
sd=vdc^0.5 
loglik=sum(log(dnorm(x=data$y, mean=pred, sd=sd)))
kc=4 
aic=-2*loglik+(2*(kc+1))
bic=-2*loglik+(log(length(res))*(kc+1))
x.plateau = knot
y.plateau = k
  tm=tmf(res,data,3)
fc=fcon2(length(res)-4,cofs,sp)
    resp = list(round(cofs, digits), round(lll, digits), knot,k,tm,fc, 
        round(r1, digits), round(r2, digits), round(aic, digits), 
        round(bic, digits), res, sres, round(shap, digits))
    names(resp) = c("Coefficients", "p values for coefficients", "x plateau","y plateau", 
        "Whole model test","Parameters, standard error and confidence intervals","R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", 
        "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")
    sin1 = ifelse(a > 0, "+", "")
    sin2 = ifelse(b > 0, "+", "")
    sin3 = ifelse(c > 0, "+", "")
    r1 = round(r1, 2)
    knot = round(knot, digits)
    e5 = substitute(atop(y == a * sin2 * b * x* sin3 * c * x, y == k) * "  " * atop(x < knot, x >= knot) * "    " * 
        R^2 * " = " * r1, list(a = a, b = b, c = c,k=k, r1 = r1, 
        knot = knot, sin2 = sin2, sin3 = sin3))
    t = list("top", "bottomright", "bottom", "bottomleft", "left", 
        "topleft", "topright", "right", "center")
    p = t[[position]]
    mei = mean
###########################################
    datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=data[,2], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)
##############################################
f1l=function(x)a+b*x+c*x^2
f2p=function(x)k+0*x
    fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
    }

ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}

bl1m=function (data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...)
{
	sdd=sd   
	names(data)=c("x","blocks","y")
	data$blocks=as.factor(data$blocks)
	mini=mean(data$x, na.rm=TRUE)
	mini=(mini/100)
	knot=seq(from=min(data$x), to=max(data$x), by=mini)
	fm=function(knot){m1=lm(y~x+(pmax(0,x-knot, na.rm=TRUE)), data=data);r=sum(resid(m1)^2);return(r)}
	s1=sapply(knot, fm)
	s=data.frame(knot,s1)
	s=s[order(s$s1),]
	knot=s[1,1]
	d=data.frame(data,pm=(pmax(0,data$x-knot, na.rm=TRUE)))
#####################
	m=lme(y~x+pm,random=~1|blocks, data=d)
	a=fixef(m)[1]
	b=fixef(m)[2]
	a=round(a,digits)
	b=round(b,digits)
	f1=function(x){a+b*x}
	d=(fixef(m)[2]+fixef(m)[3])
	resf1=f1(knot)
	c=resf1-(d*knot)
	c=round(c,digits)
	d=round(d,digits)
	f2=function(x){c+d*x}
	cofs=as.data.frame(summary(m)[20])
	cofs=cofs[,c(1,5)]
	colnames(cofs)=c("coefficientes","p_values")
	rownames(cofs)=NULL;cofs=as.data.frame(cofs)
	l1=c(a,b);l2=c(c,d)
	l=data.frame(l1,l2)
	names(l)=c("first line","second line")
	rownames(l)=c("a","b")
	l=t(l)
	k=knot
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
	r1=R2mm(m)
	r2=R3mm(m)

	res=resid(m)
	shap=shapiro.test(resid(m))
	shap=shap$p.value
	sres=scale(res)
	sres=sres[,1]
	tm=fwm(m)
	fc=fcon(m)

	resp=list(round(cofs,digits), round(l,digits), knot, tm, fc,round(r1,digits),round(r2,digits),round(AIC(m),digits),round(BIC(m),digits),res, sres, round(shap,digits) )
	names(resp)=c("Coefficients", "Lines","Knot (Break-Point)","Whole model test","Parameters, standard error and confidence intervals", "R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")

		sin1 = ifelse(a > 0, "+", "")
		sin2 = ifelse(b > 0, "+", "")
		sin3 = ifelse(c > 0, "+", "")
		sin4 = ifelse(d > 0, "+", "")

	r1=round(r1,2)
	knot=round(knot,digits)

e5 = substitute(atop(y==a*sin2*b*x,y==c*sin4*d*x)*"  "*atop(x<knot,x>=knot)*"    "* R^2 * " = " * r1, list(a = a, b = b, c = c,d = d, r1 = r1,knot=knot, sin2 = sin2, sin4 = sin4))

t = list("top", "bottomright", "bottom", "bottomleft", 
             "left", "topleft", "topright", "right", "center")
    p = t[[position]]


mei=mean
 
 datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=if(ncol(data)==2)data[,2] else data[,3], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
    
  
#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)


fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
    }
    
ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}




bl2m=function (data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...) 
{
 sdd=sd 
names(data) = c("x","blocks", "y")
data$blocks=as.factor(data$blocks)
fo=fol(data)
dl=fo[[1]]
dp=fo[[2]]
data=fo[[3]]

fml=function(data){m=lm(y~x+blocks, data=data, contrasts=list(blocks=contr.sum)); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
fmp=function(data){m=lm(y~1+blocks, data=data, contrasts=list(blocks=contr.sum)); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
om=1:length(m1)
mr=data.frame(unlist(m1), unlist(m2));names(mr)=c("l","p") 
s=rowSums(mr)
mr=data.frame(mr,s, om)
mr = mr[order(mr$s), ]
bm=mr[1,4]
bm2=mr[2,4]
fml=function(data){m=lm(y~x+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
fmp=function(data){m=lm(y~1+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m2)[[1]]
knot=(c-a)/b
bm=ifelse(knot>min(d2[,1]),bm2,bm)
fml=function(data){m=lm(y~x+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
fmp=function(data){m=lm(y~1+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
m1=lapply(dl, fml)
m2=lapply(dp, fmp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m2)[[1]]
knot=(c-a)/b
    cofs = c(a,b,c); names(cofs)=c("a","b","plateau")
    res = c(as.numeric(resid(m1)),as.numeric(resid(m2)))
    shap = shapiro.test(res)
    shap = shap$p.value
    sres = scale(res)
    sres = sres[, 1]
    vd=sum(res^2)/(length(res)-(3+(nlevels(data$blocks)-1)))

gle=(length(res)-(3+(nlevels(data$blocks)-1)))

ma1=model.matrix(m1)
maa1=cbind(ma1[,c(1,2)], c1=rep(0, length(ma1[,1])), ma1[,-c(1,2)])
ma2=model.matrix(m2)
c2=rep(0,length(ma2[,1])*2);c2=matrix(c2,ncol=2)
maa2=cbind(c2,ma2)

x=rbind(maa1,maa2)
xlin=t(x)
xx=xlin%*%x
ixx=solve(xx)
vco=ixx*vd; sp=diag(vco)^0.5
ta=a/sp[1];tb=b/sp[2]; tc=c/sp[3]
sp=c(sp[1],sp[2],sp[3])
pa=pt(abs(ta),(length(res)-(3+(nlevels(data$blocks)-1))), lower.tail = FALSE)
pa=pa*2
pb=pt(abs(tb),(length(res)-(3+(nlevels(data$blocks)-1))), lower.tail = FALSE)
pb=pb*2
pc=pt(abs(tc),(length(res)-(3+(nlevels(data$blocks)-1))), lower.tail = FALSE)
pc=pc*2
lll=c(pa,pb,pc)
names(lll)=c("a","b","plateau")
sqt = sum((data$y-mean(data$y))^2)
sqr=sqt-sum(res^2)
r1 = sqr/sqt
gl=length(res)-1
p1 =(gl/((gl + 1)-(3))*(1-r1))
r2 =1-p1
pred=na.exclude(data$y)-res
vdc=sum(res^2)/(length(res))    
sd=vdc^0.5 
loglik=sum(log(dnorm(x=data$y, mean=pred, sd=sd)))
k=3 
aic=-2*loglik+(2*(k+1))
bic=-2*loglik+(log(length(res))*(k+1-nlevels(data$blocks)))
x.plateau = knot
y.plateau = c
  tm=tm2(3, xlin,ixx,res)
fc=fcon2((length(res)-(3+(nlevels(data$blocks)-1))),cofs,sp)

    resp = list(round(cofs, digits), round(lll, digits), knot,c,tm, 
        fc,round(r1, digits), round(r2, digits), round(aic, digits), 
        round(bic, digits), res, sres, round(shap, digits))
    names(resp) = c("Coefficients", "p values for coefficients", "x plateau","y plateau", 
        "Whole model test","Parameters, standard error and confidence intervals", "R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", 
        "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")
    sin1 = ifelse(a > 0, "+", "")
    sin2 = ifelse(b > 0, "+", "")
    sin3 = ifelse(c > 0, "+", "")
    r1 = round(r1, 2)
    knot = round(knot, digits)
    e5 = substitute(atop(y == a * sin2 * b * x, y == c) * "  " * atop(x < knot, x >= knot) * "    " * 
        R^2 * " = " * r1, list(a = a, b = b, c = c, r1 = r1, 
        knot = knot, sin2 = sin2))
    t = list("top", "bottomright", "bottom", "bottomleft", "left", 
        "topleft", "topright", "right", "center")
    p = t[[position]]
    mei = mean


###########################################
    
    datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=if(ncol(data)==2)data[,2] else data[,3], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
    

#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)
##############################################

f1l=function(x)a+b*x
f2p=function(x)c+0*x
    fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl,xlab = xlab,ylab=ylab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }

    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao){
         plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty)
    }
    
ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}


bl3m=function(data, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,...)
{
sdd=sd       
names(data) = c("x", "blocks","y")
data$blocks=as.factor(data$blocks)
fo=fop(data)
dl=fo[[1]]
dp=fo[[2]]
data=fo[[3]]
   
fmq=function(data){m=lm(y~x+I(x^2)+blocks, data=data, contrasts=list(blocks=contr.sum));r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
fmp=function(data){m=lm(y~1+blocks, data=data,contrasts=list(blocks=contr.sum)); r1=as.numeric(resid(m)); r=sum(r1^2);return(r)}
m1=lapply(dl, fmq)
m2=lapply(dp, fmp)
om=1:length(m1)
mr=data.frame(unlist(m1), unlist(m2));names(mr)=c("l","p") 
s=rowSums(mr)
mr=data.frame(mr,s, om)
mr = mr[order(mr$s), ]
bm=mr[1,4]
mq=function(data){m=lm(y~x+I(x^2)+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
m1=lapply(dl, mq)
mp=function(data){m=lm(y~1+blocks, data=data, contrasts=list(blocks=contr.sum)); return(m)}
m2=lapply(dp, mp)
m1=m1[[bm]]
m2=m2[[bm]]
datao=data
d1=dl[[bm]]
d2=dp[[bm]]
a=coef(m1)[[1]]
b=coef(m1)[[2]]
c=coef(m1)[[3]]
k=coef(m2)[[1]]
delta=sqrt((b^2)-(4*c*a)+(4*c*k))
knot1=((-1*b)+delta)/(2*c)
knot2=((-1*b)-delta)/(2*c)
k1=abs(knot1)
k2=abs(knot2)
kk=c(k1,k2); knot=sort(kk)[1]
    cofs = c(a,b,c,k); names(cofs)=c("a","b","c","plateau")
    res = c(as.numeric(resid(m1)[1:length(d1$x)]),as.numeric(resid(m2)))
    shap = shapiro.test(res)
    shap = shap$p.value
    sres = scale(res)
    sres = sres[, 1]
    glb=nlevels(data$blocks)-1
    vd=sum(res^2)/(length(res)-4-glb)

ma1=model.matrix(m1)
maa1=cbind(ma1[,c(1,2,3)], c1=rep(0, length(ma1[,1])), ma1[,-c(1,2,3)])
ma2=model.matrix(m2)
c2=rep(0,length(ma2[,1])*3);c2=matrix(c2,ncol=3)
maa2=cbind(c2,ma2)

x=rbind(maa1,maa2)
xlin=t(x)
xx=xlin%*%x
ixx=solve(xx)
vco=ixx*vd; sp=diag(vco)^0.5
ta=a/sp[1];tb=b/sp[2]; tc=c/sp[3]; tk=k/sp[4]
pa=pt(abs(ta),(length(res)-4-glb), lower.tail = FALSE)
pa=pa*2
pb=pt(abs(tb),(length(res)-4-glb), lower.tail = FALSE)
pb=pb*2
pc=pt(abs(tc),(length(res)-4-glb), lower.tail = FALSE)
pc=pc*2
pk=pt(abs(tk),(length(res)-4-glb), lower.tail = FALSE)
pk=pk*2
lll=c(pa,pb,pc, pk)
names(lll)=c("a","b","c","plateau")
sqt = sum((data$y-mean(data$y))^2)
sqr=sqt-sum(res^2)
r1 = sqr/sqt
gl=length(res)-1
p1 =(gl/((gl + 1)-(4))*(1-r1))
r2 =1-p1
pred=data$y-res
vdc=sum(res^2)/(length(res))    
sd=vdc^0.5 
loglik=sum(log(dnorm(x=data$y, mean=pred, sd=sd)))
kc=4
aic=-2*loglik+(2*(kc+1))
bic=-2*loglik+(log(length(res)))*(kc+1-nlevels(data$blocks))
x.plateau = knot
y.plateau = k
  tm=tm2(4, xlin,ixx,res)
fc=fcon2(length(res)-4,cofs,sp[c(1,2,3,4)])
    resp = list(round(cofs, digits), round(lll, digits), knot,k,tm,fc, 
        round(r1, digits), round(r2, digits), round(aic, digits), 
        round(bic, digits), res, sres, round(shap, digits))
    names(resp) = c("Coefficients", "p values for coefficients", "x plateau","y plateau", 
        "Whole model test","Parameters, standard error and confidence intervals","R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", 
        "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")
    sin1 = ifelse(a > 0, "+", "")
    sin2 = ifelse(b > 0, "+", "")
    sin3 = ifelse(c > 0, "+", "")
    r1 = round(r1, 2)
    knot = round(knot, digits)
    e5 = substitute(atop(y == a * sin2 * b * x* sin3 * c * x, y == k) * "  " * atop(x < knot, x >= knot) * "    " * 
        R^2 * " = " * r1, list(a = a, b = b, c = c,k=k, r1 = r1, 
        knot = knot, sin2 = sin2, sin3 = sin3))
    t = list("top", "bottomright", "bottom", "bottomleft", "left", 
        "topleft", "topright", "right", "center")
    p = t[[position]]
    mei = mean
###########################################
   datao = means(data)
    dataoi = ifelse(mei == TRUE, 1, 2)
    data=data.frame(x=data[,1], y=if(ncol(data)==2)data[,2] else data[,3], sd=rep(0,length(data[,1])))
    ldi = list(datao, data)
    datao = ldi[[dataoi]]
#################################################
    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)
##############################################
f1l=function(x)a+b*x+c*x^2
f2p=function(x)k+0*x
    fp1 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp2 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }
fp3 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
    fp4 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot,  col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
             col = col, lty = lty,...)
    }

fp11 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
        legend(p, legend = e5, bty = "n", cex = 0.8)
    }
fp22 = function(datao) {
        plot(datao[, 2] ~ as.numeric(as.character(datao[, 1])), 
            data = datao, ylim = yl, xlim = xl, ylab = ylab, xlab = xlab, bty = "n", pch=pch,...)
        curve(f1l, add = TRUE, from = min(data$x, na.rm = TRUE), 
            to = knot, col = col, lty = lty,...)
        curve(f2p, add = TRUE, from = knot, to = max(data$x, na.rm = TRUE), 
            col = col, lty = lty,...)
    }

ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
    ll[[li]](datao)
    return(resp)
}



ggg=list(bl1,bl2,bl3,bl1m,bl2m,bl3m)
ggg[[model]](data, xlab = xlab, ylab = ylab, 
    position = position, digits = digits, mean = mean, sd=sd, legend = legend, lty=lty, col=col, pch=pch,...)

}




