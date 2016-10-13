bl=function(data, xlab="Explanatory Variable", ylab="Response Variable", position=1, digits=6, mean=TRUE, legend=TRUE)
{

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
resp=list(round(cofs,digits), round(l,digits), knot, round(r1,digits),round(r2,digits),round(AIC(m),digits),round(BIC(m),digits),res, sres, round(shap,digits) )
names(resp)=c("Coefficients", "Lines","Knot (Break-Point)", "R-squared", "Adjusted R-squared", "AIC", "BIC", "Residuals", "Standartized residuals", "P value (Shapiro-Wilk test for residuals)")

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
means = function(data) {
        t = as.factor(data[, 1])
        d = data.frame(t, data[, -1])
        s = split(data.frame(d[, -1]), d$t)
        r = lapply(s, colMeans, na.rm = TRUE)
        r = lapply(r, round, 2)
        rr = t(data.frame(r))
        rr = data.frame(rr)
        rownames(rr) = NULL
        treat = levels(t)
        rr = data.frame(treat, rr)
        colnames(rr) = colnames(data)
        return(rr)
    }

datao=means(data)
dataoi=ifelse(mei==TRUE,1,2)
	ldi=list(datao,data)
	datao=ldi[[dataoi]]


	minx = min(data[, 1], na.rm=TRUE) - sd(data[, 1], na.rm=TRUE)/2
	maxx = max(data[, 1], na.rm=TRUE) + sd(data[, 1], na.rm=TRUE)/2
        miny = min(data[, 2], na.rm=TRUE) - sd(data[, 2], na.rm=TRUE)/2
        maxy = max(data[, 2], na.rm=TRUE) + sd(data[, 2], na.rm=TRUE)/2


li=ifelse(legend==TRUE,1,2)


fp1=function(datao){
plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, ylim=c(miny,maxy), xlim=c(minx,maxx), ylab=ylab, xlab=xlab, bty="n")
curve(f1,add=TRUE, from=min(data$x, na.rm=TRUE), to=knot, col = "dark red", lty = 2)
curve(f2, add=TRUE, from=knot, to=max(data$x, na.rm=TRUE), col = "dark blue", lty = 3)
legend(p,legend=e5,bty = "n", cex=0.8)
}

fp2=function(datao){
plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, ylim=c(miny,maxy), xlim=c(minx,maxx), ylab=ylab, xlab=xlab, bty="n")
curve(f1,add=TRUE, from=min(data$x, na.rm=TRUE), to=knot, col = "dark red", lty = 2)
curve(f2, add=TRUE, from=knot, to=max(data$x, na.rm=TRUE), col = "dark blue", lty = 3)
}

ll=list(fp1,fp2)

ll[[li]](datao)

return(resp)

}
