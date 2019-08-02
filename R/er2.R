er2<-
function(data, design=1, list=FALSE, type=2){
    list=ifelse(list==FALSE, 1,2)
    type=ifelse(type==1, 1,2)
fa11=function(a){
        a=anova(a); res=a;d=data.frame(res); d=round(d,4); d1=d[,5]; d2=ifelse(d1<0.001, "<0.001", d1); 
        d2=d2[-length(d2)];d2=c(d2,"-"); d=d[,-5];d=data.frame(d,d2);d[is.na(d)] <- "-"
        names(d)=c("df", "type I SS", "mean square", "F value", "p>F")
        return(d)
    }

 fa22=function(m){
	a=drop1(m,.~.,test="F"); ress=c(m$df.residual,sum(m$residuals^2),sum(m$residuals^2)/m$df.residual,NA,NA)
	res=a; d=data.frame(res[-1,]); d=data.frame(d[,1],d[,2],d[,2]/d[,1],d[,5],d[,6])
        ;d=rbind(d,ress);d=round(d,4);d1=d[,5]; d2=ifelse(d1<0.001, "<0.001", d1);
        d2=d2[-length(d2)];d2=c(d2,"-"); d=d[,-5];d=data.frame(d,d2);d[is.na(d)] <- "-"
        names(d)=c("df", "type III SS", "mean square", "F value", "p>F"); rownames(d)=c(rownames(res[-1,]),"residuals")
        return(d)
    }


lfun=list(fa11,fa22)
fa1=lfun[[type]]
            f1=function(data){ 
                names(data)=c("treatments","response")
                data=data.frame(treatments=as.numeric(data$treatments), response=data$response, linear=data$treatments, quadratic=data$treatments^2, cubic=data$treatments^3, treatment=factor(data$treatments) )
                m1<-lm(response~linear,data=data)
                m2<-lm(response~linear+quadratic,data=data)
                m3<-lm(response~linear+quadratic+cubic,data=data)
                anova=list(fa1(m1),fa1(m2),fa1(m3))
                names(anova)=c('linear','quadratic','cubic')
                m11=lm(response~linear, data=data)
                m12=lm(response~linear+quadratic, data=data)
                m13=lm(response~linear+quadratic+cubic, data=data)

                models=list(round(m11[1][[1]],4), round(m12[1][[1]],4), round(m13[1][[1]],4))
                names(models)=c('linear','quadratic','cubic')
                c1=round(summary(m1)[[4]],4)
                c2=round(summary(m2)[[4]],4)
                c3=round(summary(m3)[[4]],4)
                c1=as.data.frame(c1);c1=c1[c(1,2),]
                c2=as.data.frame(c2);c2=c2[c(1,2,3),]
                c3=as.data.frame(c3);c3=c3[c(1,2,3,4),]
                
                coefs=list(c1,c2,c3)
                names(coefs)=c('linear','quadratic','cubic')
                st=sum(anova(m1)[[2]][c(1,2)])
                st1=sum(anova(m1)[[2]][1])
                st2=sum(anova(m2)[[2]][c(1,2)])
                st3=sum(anova(m3)[[2]][c(1,2,3)])
                r1=summary(m1)$r.squared*100
                r2=summary(m2)$r.squared*100
                r3=summary(m3)$r.squared*100
                r12=summary(m1)$adj.r.squared*100
                r22=summary(m2)$adj.r.squared*100
                r23=summary(m3)$adj.r.squared*100
		AIC=round(c(AIC(m1),AIC(m2),AIC(m3)),4)
		BIC=round(c(BIC(m1),BIC(m2),BIC(m3)),4)
                R_squared=round(c(r1,r2,r3),4)
                Adjusted_R_squared=round(c(r12,r22,r23),4)
                Models=c('linear','quadratic','cubic')
                r=data.frame(Models,R_squared, Adjusted_R_squared, AIC, BIC)
                rf1=list(anova,models,coefs,r)
                rf2= list(anova[c(1,2)],models[c(1,2)],coefs[c(1,2)],r[c(1,2),])
                rf=list(rf1,rf2)
                j=nlevels(data$treatment)
                j=ifelse(j>3,1,2)
                rf=rf[[j]]
                names(rf)=c("Analysis of variance","Models","t test for coefficients","R-squared")
                return(rf)
            }
            
            f2=function(data){ 
                names(data)=c("treatments","blocks","response")
                data=data.frame(treatments=as.numeric(data$treatments), blocks=as.factor(data$blocks), response=data$response, linear=data$treatments, quadratic=data$treatments^2, cubic=data$treatments^3, treatment=factor(data$treatments) )
                m1<-lm(response~linear+blocks,data=data, contrasts=list(blocks=contr.sum))
                m2<-lm(response~linear+quadratic+blocks,data=data, contrasts=list(blocks=contr.sum))
                m3<-lm(response~linear+quadratic+cubic+blocks,data=data, contrasts=list(blocks=contr.sum))
                anova=list(fa1(m1),fa1(m2),fa1(m3))
                names(anova)=c('linear','quadratic','cubic')
                m11=lm(response~linear, data=data)
                m12=lm(response~linear+quadratic, data=data)
                m13=lm(response~linear+quadratic+cubic, data=data)
                
                models=list(round(m11[1][[1]],4), round(m12[1][[1]],4), round(m13[1][[1]],4))
                names(models)=c('linear','quadratic','cubic')
                c1=round(summary(m1)[[4]],4)
                c2=round(summary(m2)[[4]],4)
                c3=round(summary(m3)[[4]],4)
                c1=as.data.frame(c1);c1=c1[c(1,2),]
                c2=as.data.frame(c2);c2=c2[c(1,2,3),]
                c3=as.data.frame(c3);c3=c3[c(1,2,3,4),]
                coefs=list(c1,c2,c3)
                names(coefs)=c('linear','quadratic','cubic')
                st=sum(anova(m1)[[2]][c(1,2)])
                st1=sum(anova(m1)[[2]][1])
                st2=sum(anova(m2)[[2]][c(1,2)])
                st3=sum(anova(m3)[[2]][c(1,2,3)])
                r1=summary(m1)$r.squared*100
                r2=summary(m2)$r.squared*100
                r3=summary(m3)$r.squared*100
                r12=summary(m1)$adj.r.squared*100
                r22=summary(m2)$adj.r.squared*100
                r23=summary(m3)$adj.r.squared*100
		AIC=round(c(AIC(m1),AIC(m2),AIC(m3)),4)
		BIC=round(c(BIC(m1),BIC(m2),BIC(m3)),4)
                R_squared=round(c(r1,r2,r3),4)
                Adjusted_R_squared=round(c(r12,r22,r23),4)
                Models=c('linear','quadratic','cubic')
                r=data.frame(Models,R_squared, Adjusted_R_squared, AIC, BIC)
                rf1=list(anova,models,coefs,r)
                rf2= list(anova[c(1,2)],models[c(1,2)],coefs[c(1,2)],r[c(1,2),])
                rf=list(rf1,rf2)
                j=nlevels(data$treatment)
                j=ifelse(j>3,1,2)
                rf=rf[[j]]
                names(rf)=c("Analysis of variance","Models","t test for coefficients","R-squared")
                return(rf)
            }
             f3=function(data){ 
                names(data)=c("treatments","rows", "columns","response")
                data=data.frame(treatments=as.numeric(data$treatments), rows=as.factor(data$rows), columns=as.factor(data$columns),  response=data$response, linear=data$treatments, quadratic=data$treatments^2, cubic=data$treatments^3, treatment=factor(data$treatments) )
                         m1<-lm(response~linear+rows+columns,data=data, contrasts=list(rows=contr.sum, columns=contr.sum))
                                m2<-lm(response~linear+quadratic+rows+columns,data=data, contrasts=list(rows=contr.sum, columns=contr.sum))
                                  m3<-lm(response~linear+quadratic+cubic+rows+columns,data=data, contrasts=list(rows=contr.sum, columns=contr.sum))
                     anova=list(fa1(m1),fa1(m2),fa1(m3))
                names(anova)=c('linear','quadratic','cubic')
                m11=lm(response~linear, data=data)
                m12=lm(response~linear+quadratic, data=data)
                m13=lm(response~linear+quadratic+cubic, data=data)
                models=list(round(m11[1][[1]],4), round(m12[1][[1]],4), round(m13[1][[1]],4))
                names(models)=c('linear','quadratic','cubic')
                c1=round(summary(m1)[[4]],4)
                c2=round(summary(m2)[[4]],4)
                c3=round(summary(m3)[[4]],4)
                c1=as.data.frame(c1);c1=c1[c(1,2),]
                c2=as.data.frame(c2);c2=c2[c(1,2,3),]
                c3=as.data.frame(c3);c3=c3[c(1,2,3,4),]
                coefs=list(c1,c2,c3)
                names(coefs)=c('linear','quadratic','cubic')
                st=sum(anova(m1)[[2]][c(1,2)])
                st1=sum(anova(m1)[[2]][1])
                st2=sum(anova(m2)[[2]][c(1,2)])
                st3=sum(anova(m3)[[2]][c(1,2,3)])
                r1=summary(m1)$r.squared*100
                r2=summary(m2)$r.squared*100
                r3=summary(m3)$r.squared*100
                r12=summary(m1)$adj.r.squared*100
                r22=summary(m2)$adj.r.squared*100
                r23=summary(m3)$adj.r.squared*100
		AIC=round(c(AIC(m1),AIC(m2),AIC(m3)),4)
		BIC=round(c(BIC(m1),BIC(m2),BIC(m3)),4)
                R_squared=round(c(r1,r2,r3),4)
                Adjusted_R_squared=round(c(r12,r22,r23),4)
                Models=c('linear','quadratic','cubic')
                r=data.frame(Models,R_squared, Adjusted_R_squared, AIC, BIC)
                rf1=list(anova,models,coefs,r)
                rf2= list(anova[c(1,2)],models[c(1,2)],coefs[c(1,2)],r[c(1,2),])
                rf=list(rf1,rf2)
                j=nlevels(data$treatment)
                j=ifelse(j>3,1,2)
                rf=rf[[j]]
                names(rf)=c("Analysis of variance","Models","t test for coefficients","R-squared")
                return(rf)
            }
           
            f4=function(data){ 
                names(data)=c("treatments", "squares","rows", "columns","response")
                data=data.frame(treatments=as.numeric(data$treatments), rows=as.factor(data$rows), squares=as.factor(data$squares), columns=as.factor(data$columns),  response=data$response, linear=data$treatments, quadratic=data$treatments^2, cubic=data$treatments^3, treatment=factor(data$treatments) )
                    m1<-lm(response~linear+squares+rows+columns,data=data, contrasts=list(squares=contr.sum, rows=contr.sum, columns=contr.sum))
                    m2<-lm(response~linear+quadratic+squares+rows+columns,data=data, contrasts=list(squares=contr.sum,rows=contr.sum, columns=contr.sum))
                              m3<-lm(response~linear+quadratic+cubic+squares+rows+columns,data=data, contrasts=list(squares=contr.sum,rows=contr.sum, columns=contr.sum))
                                    
                anova=list(fa1(m1),fa1(m2),fa1(m3))
                names(anova)=c('linear','quadratic','cubic')
                m11=lm(response~linear, data=data)
                m12=lm(response~linear+quadratic, data=data)
                m13=lm(response~linear+quadratic+cubic, data=data)
                models=list(round(m11[1][[1]],4), round(m12[1][[1]],4), round(m13[1][[1]],4))
                names(models)=c('linear','quadratic','cubic')
                c1=round(summary(m1)[[4]],4)
                c2=round(summary(m2)[[4]],4)
                c3=round(summary(m3)[[4]],4)
                c1=as.data.frame(c1);c1=c1[c(1,2),]
                c2=as.data.frame(c2);c2=c2[c(1,2,3),]
                c3=as.data.frame(c3);c3=c3[c(1,2,3,4),]
                coefs=list(c1,c2,c3)
                names(coefs)=c('linear','quadratic','cubic')
                st=sum(anova(m1)[[2]][c(1,2)])
                st1=sum(anova(m1)[[2]][1])
                st2=sum(anova(m2)[[2]][c(1,2)])
                st3=sum(anova(m3)[[2]][c(1,2,3)])
                r1=summary(m1)$r.squared*100
                r2=summary(m2)$r.squared*100
                r3=summary(m3)$r.squared*100
                r12=summary(m1)$adj.r.squared*100
                r22=summary(m2)$adj.r.squared*100
                r23=summary(m3)$adj.r.squared*100
		AIC=round(c(AIC(m1),AIC(m2),AIC(m3)),4)
		BIC=round(c(BIC(m1),BIC(m2),BIC(m3)),4)
                R_squared=round(c(r1,r2,r3),4)
                Adjusted_R_squared=round(c(r12,r22,r23),4)
                Models=c('linear','quadratic','cubic')
                r=data.frame(Models,R_squared, Adjusted_R_squared, AIC, BIC)
                rf1=list(anova,models,coefs,r)
                rf2= list(anova[c(1,2)],models[c(1,2)],coefs[c(1,2)],r[c(1,2),])
                rf=list(rf1,rf2)
                j=nlevels(data$treatment)
                j=ifelse(j>3,1,2)
                rf=rf[[j]]
                names(rf)=c("Analysis of variance","Models","t test for coefficients","R-squared")
                return(rf)
            }
      
de1 = c(1)
    de2 = c(1, 2)
    de3 = c(1, 2, 3)
    de4 = c(1, 2, 3, 4)
de = list(de1, de2, de3, de4)
 de = de[[design]]
    d = as.list(data)
    d1 = d[de]
    d2 = d[-de]
    f = function(h) {
        data.frame(d1, d2[h])
    }
    h = length(d2)
    h = 1:h
    l = lapply(h, f)
    l2 = list(f1, f2, f3, f4)
    fun = l2[[design]]
    li1 = lapply(l, fun)
    names(li1) = names(d2)
    li = list(fun(data), li1)
    li = li[[list]]
    return(li)
}
