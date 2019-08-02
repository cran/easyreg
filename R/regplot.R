
regplot <-
function(data, model=1, start=c(a=1,b=1,c=1,d=1,e=1), xlab="Explanatory Variable", ylab="Response Variable", position=1, digits=6, mean=TRUE, sd=FALSE, legend = TRUE, lty=2, col="dark blue", pch=20,xlim="defalt.x",ylim="defalt.y",...)

{
mei=mean
li=legend
lty=lty
col=col
pch=pch
sdd=sd

	minx = min(data[, 1], na.rm=TRUE) - sd(data[, 1],na.rm=TRUE)/2
        maxx = max(data[, 1],na.rm=TRUE) + sd(data[, 1],na.rm=TRUE)/2
        miny = min(data[, 2],na.rm=TRUE) - sd(data[, 2],na.rm=TRUE)/2
        maxy = max(data[, 2],na.rm=TRUE) + sd(data[, 2],na.rm=TRUE)/2
      
t1=min(data[,2], na.rm=TRUE)
t2=max(data[,2], na.rm=TRUE)  

xli=ifelse(xlim=="defalt.x",1,2)
xli=xli[1]
xll=list(c(minx,maxx),xlim)
xl=xll[[xli]]

yli=ifelse(ylim=="defalt.y",1,2)
yli=yli[1]
yll=list(c(miny,maxy),ylim)
yl=yll[[yli]]



    rest=er1(data=data, model=model, start=c(a=start[1],b=start[2],c=start[3],d=start[4],e=start[5]), digits=digits)
    
    res=rest[[1]][[1]]

    
    
    t = list("top", "bottomright", "bottom", "bottomleft", 
             "left", "topleft", "topright", "right", "center")
    p = t[[position]]
    
        
iii=1:(ncol(data)-1)
fr=function(iii){
res2=rest[[iii]][[1]];return(res2)}

lr=lapply(iii,fr)
res=data.frame(lr);names(res)=colnames(data)[-1]

 plot1=function(data)
    {
        f1=function(i)res[1,i]+res[2,i]*se
        f2=function(i)res[1,i]+res[2,i]*se+res[3,i]*se^2
        f3=function(i)res[1,i] + res[2,i] * (se - res[3,i]) * (se <= res[3,i])
        f4=function(i)(res[1,i] + res[2,i] * se + res[3,i] * I(se^2)) * (se <= -0.5 * res[2,i]/res[3,i]) + 
            (res[1,i] + I(-res[2,i]^2/(4 * res[3,i]))) * (se > -0.5 * res[2,i]/res[3,i])
        f5=function(i)ifelse(se>=res[4,i],(res[1,i]-res[3,i]*res[4,i])+(res[2,i] +res[3,i])*se, res[1,i]+res[2,i] *se)
        f6=function(i)res[1,i]*exp(res[2,i]*se)
        f7=function(i)res[1,i]*(1+res[2,i]*(exp(-res[3,i]*se)))^-1
        f8=function(i)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))^3
        f9=function(i)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))
        f10=function(i)res[1,i]*exp(-res[2,i]*exp(-res[3,i]*se))
        f11=function(i)(res[1,i]*se^res[2,i])*exp(-res[3,i]*se)
        f12=function(i)res[1,i] + res[2,i] * (1 - exp(-res[3,i] * se))
        f13=function(i)(res[1,i]/(1+exp(2-4*res[3,i]*(se-res[5,i]))))+(res[2,i]/(1+exp(2-4*res[4,i]*(se-res[5,i]))))
        f14=function(i)res[1,i]*(se^res[2,i])
 	f15=function(i)res[1,i]+res[2,i]*se+res[3,i]*se^2+res[4,i]*se^3
	f16=function(i)res[1,i]/(1+res[2,i]*(exp(-res[3,i]*se)))^res[4,i]
	f17=function(i)(res[1,i]^res[4,i]+ ((res[2,i]^res[4,i])-(res[1,i]^res[4,i]) )*((1-exp(-res[3,i]*(se-t1)))/ (1-exp(-res[3,i]*(t2-t1)))))^(1/res[4,i])
        mod=list(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17)
        se=seq(min(data[,1]),max(data[,1]),by=0.01)
        i=1:(ncol(data)-1)
        pred=lapply(i,mod[[model]])
        u=unlist(pred)
        mat=matrix(u, ncol=length(pred))
        pred=as.data.frame(mat)
        dd=data.frame(se,pred)
        n=names(data)
        names=n[-1]
        matplot(dd[,1], dd[,-1], type = "l", bty = "n", xlab = xlab, 
                ylab = ylab, col = c(1:(ncol(data)-1)), lty = c(1:(ncol(data)-1)))
        legend(p, names, bty = "n", col = c(1:(ncol(data)-1)), 
               lty = c(1:(ncol(data)-1)))
    }
    



    plot2=function(data)
    {
        ff1=function(se)res[1,i]+res[2,i]*se
        ff2=function(se)res[1,i]+res[2,i]*se+res[3,i]*se^2
        ff3=function(se)res[1,i] + res[2,i] * (se - res[3,i]) * (se <= res[3,i])
        ff4=function(se)(res[1,i] + res[2,i] * se + res[3,i] * I(se^2)) * (se <= -0.5 * res[2,i]/res[3,i]) + 
            (res[1,i] + I(-res[2,i]^2/(4 * res[3,i]))) * (se > -0.5 * res[2,i]/res[3,i])
        ff5=function(se)ifelse(se>=res[4,i],(res[1,i]-res[3,i]*res[4,i])+(res[2,i] +res[3,i])*se, res[1,i]+res[2,i] *se)
        ff6=function(se)res[1,i]*exp(res[2,i]*se)
        ff7=function(se)res[1,i]*(1+res[2,i]*(exp(-res[3,i]*se)))^-1
        ff8=function(se)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))^3
        ff9=function(se)res[1,i]*(1-res[2,i]*(exp(-res[3,i]*se)))
        ff10=function(se)res[1,i]*exp(-res[2,i]*exp(-res[3,i]*se))
        ff11=function(se)(res[1,i]*se^res[2,i])*exp(-res[3,i]*se)
        ff12=function(se)res[1,i] + res[2,i] * (1 - exp(-res[3,i] * se))
        ff13=function(se)  (res[1,i]/(1+exp(2-4*res[3,i]*(se-res[5,i]))))+(res[2,i]/(1+exp(2-4*res[4,i]*(se-res[5,i]))))
        ff14=function(se)res[1,i]*(se^res[2,i])
 	ff15=function(se)res[1,i]+res[2,i]*se+res[3,i]*se^2+res[4,i]*se^3
	ff16=function(se)res[1,i]/(1+res[2,i]*(exp(-res[3,i]*se)))^res[4,i]
	ff17=function(se)(res[1,i]^res[4,i]+ ((res[2,i]^res[4,i])-(res[1,i]^res[4,i]) )*((1-exp(-res[3,i]*(se-t1)))/ (1-exp(-res[3,i]*(t2-t1)))))^(1/res[4,i])
        mod2=list(ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10,ff11,ff12,ff13,ff14,ff15,ff16,ff17)
        i=1:(ncol(data)-1)
        

        c1=res[1,];c2=res[2,];c3=res[3,];c4=res[4,];c5=res[5,]

	cc1=c1-c3*c4; cc1=round(cc1,digits)
	cc2=c2+c3; cc2=round(cc2,digits)
	pmm = (c1 + I(-c2^2/(4 * c3))); pmm=round(pmm, digits)
        pcc = -0.5 * c2/c3; pcc=round(pcc, digits)
        
        sin1 = ifelse(c1 > 0, "+", "")
        sin2 = ifelse(c2 > 0, "+", "")
        sin3 = ifelse(c3 > 0, "+", "")
        sin4 = ifelse(c4 > 0, "+", "")
        sin5 = ifelse(c5 > 0, "+", "")
        
        r2=res["r-squared",]
        
        e1 = substitute(y == c1 * sin2 * c2 * x * "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                r2 = r2, sin2 = sin2))
        
        e2 = substitute(y == c1 * sin2 * c2 * x * sin3 * c3 * 
                            x^2 * "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                c3 = c3, r2 = r2, sin2 = sin2, sin3 = sin3))
        e3 = substitute(atop(y == c1 * sin2 * c2 * (x -c3),y==c1*"  plateau")*"  "*atop(x<c3,x>=c3) * 
                            "  " * R^2 *" = " * r2, list(c1 = c1, c2 = c2, r2 = r2, c3 = c3, 
                                                         sin2 = sin2))
        e4 = substitute(atop(y == c1 * sin2 * c2 * x * sin3 * c3 * 
                            x^2,y==pmm*"   plateau")* "   "*atop(x<pcc,x>=pcc) *  "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, c3 = c3, 
                                                                 r2 = r2, sin2 = sin2, sin3 = sin3, pmm=pmm, pcc=pcc))
        e5 = substitute(atop(y==c1*sin2*c2*x,y==cc1+cc2*x)*"   "*atop(x<c4,x>=c4)*"    "* R^2 * " = " * r2, list(c1 = c1, c2 = c2, cc1 = cc1,cc2 = cc2,c4=c4,r2 = r2, sin2 = sin2))
        
        e6 = substitute(y==c1*e^{c2*x}*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                     r2 = r2)) 
        e7 = substitute(y==c1*(1*sin2*c2*e^{-c3*x})^-1*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                     sin2 = sin2, c3 = c3, r2 = r2)) 
        e8 = substitute(y==c1*(1-c2*e^{-c3*x})^3*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                               c3 = c3, r2 = r2)) 
        e9 = substitute(y==c1*(1-c2*e^{-c3*x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                             c3 = c3, r2 = r2)) 
        e10 = substitute(y==c1*e^(-c2*e^{-c3*x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                               c3 = c3, r2 = r2)) 
        e11 = substitute(y==c1*x^c2*e^{-c3*x}*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                            c3 = c3, r2 = r2))
        e12 = substitute(y==c1 + c2 * (1 - e^{-c3 * x})*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                      c3 = c3, r2 = r2))
        cc3=4*c3
        cc4=4*c4
        e13 = substitute(y==frac(c1,1+e^(2-cc3*(x-c5)))+frac(c2,1+e^(2-cc4*(x-c5)))*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                                                  cc3 = cc3, r2 = r2, cc4=cc4, c5=c5))
 	e14 = substitute(y==c1*{x^c2}*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                     r2 = r2)) 
        e15 = substitute(y == c1 * sin2 * c2 * x * sin3 * c3 * 
                            x^2 *sin4 * c4 *  "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                c3 = c3,c4=c4, r2 = r2, sin2 = sin2, sin3 = sin3, sin4=sin4))

 	e16 = substitute(y==frac(c1,(1*sin2*c2*e^{-c3*x})^c4)*"  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2, 
                                                                                     sin2 = sin2, c3 = c3,c4=c4, r2 = r2))

	e17 = substitute(y==c1^c4+(c2^c4-c1^c4)*bgroup("(",frac(1-e^-c3*(t-t1), 1-e^-c3*(t2-t1)),")")^frac(1,c4)*  "  " * R^2 * " = " * r2, list(c1 = c1, c2 = c2,c3 = c3,c4=c4, r2 = r2))
                                                               
         # demo(plotmath)
        
        ee=list(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14, e15,e16, e17)
        
        eee=ee[[model]]

    means = function(data) {
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

	

datao=means(data)
	dataoi=ifelse(mei==TRUE,1,2)
	ldi=list(datao,data)
	datao=ldi[[dataoi]]
        
fp1=function(datao){
        plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, bty="n",pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        legend(p,legend=eee,bty = "n", cex=0.8)
}

fp2=function(datao){
plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, bty="n",pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
arrows(as.numeric(as.character(datao[, 1])),datao[, 2]-datao[, 3],as.numeric(as.character(datao[, 1])),datao[, 2]+datao[, 3], code=3,length=0.04, angle=90)
        }
fp3=function(datao){
        plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, bty="n",pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
        legend(p,legend=eee,bty = "n", cex=0.8)
}

fp4=function(datao){
plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, bty="n", pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
        }
fp11=function(datao){
        plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab,  bty="n",pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
        legend(p,legend=eee,bty = "n", cex=0.8)
}

fp22=function(datao){
plot(datao[,2]~as.numeric(as.character(datao[,1])), data=datao, xlab = xlab, 
             ylab = ylab, bty="n",pch=pch,xlim=xl,ylim=yl,...)
        plot(mod2[[model]], min(data[,1]), max(data[,1]), add = TRUE, 
             lty=lty, col=col,...)
        }


    li = ifelse(legend == TRUE, 1, 2)
    lisd=ifelse(sdd==TRUE,0,2)
    li=sum(li+lisd)

ll = list(fp1, fp2, fp3, fp4, fp11, fp22)
li=ifelse(mean==TRUE, li, 5)
le=ifelse(legend==TRUE, 1,2)
lie=li+le
li=ifelse(lie==7,6,li)
ll[[li]](datao)


return(rest)

    }
 ppp=ifelse(ncol(data)==2,2,1)
    lll=list(plot1,plot2)
    lll[[ppp]](data)   
    
}
