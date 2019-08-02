\name{er1}
\alias{er1}

\title{
Analysis of regression
}
\description{
The function performs analysis of some linear and nonlinear models
}
\usage{
er1(data, model = 1, start = c(a = 1, b = 1, c = 1, d = 1, e = 1), 
mixed=FALSE, digits=6, alpha=0.05)
}

\arguments{
  \item{data}{
data is a data.frame

The first column should contain the treatments (explanatory variable) and the remaining columns the response variables.
}
  \item{model}{
define the model

1 = "y~a+b*x" linear

2 = "y~a+b*x+c*x^2" quadratic

3 = "y ~ a + b * (x - c) * (x <= c)" linear plateau

4 = "y ~ (a + b * x + c * I(x^2)) * (x <= -0.5 * b/c) + (a + I(-b^2/(4 * c))) * (x > -0.5 * b/c)" quadratic plateau
	
5 = "ifelse(x>=d,(a-c*d)+(b+c)*x, a+b*x)" two linear

6 = "y~a*exp(b*x)" exponential

7 = "y~a*(1+b*(exp(-c*x)))^-1" logistic

8 = "y~a*(1-b*(exp(-c*x)))^3" van bertalanffy

9 = "y~a*(1-b*(exp(-c*x)))" brody

10 = "y~a*exp(-b*exp(-c*x)" gompertz

11 = "y~(a*x^b)*exp(-c*x)" lactation curve

12 = "y ~ a + b * (1 - exp(-c * x))" ruminal degradation curve

13 = "y~(a/(1+exp(2-4*c*(x-e))))+(b/(1+exp(2-4*d*(x-e))))" logistic bi-compartmental

14 = "y~a*(x^b)" exponential (allometric model)

15 = "y~a+b*x+c*x^2+d*x^3" cubic

16 = "y~a/(1+b*(exp(-c*x)))^d" richards

17 = "y~(a^d+ ((b^d)-(a^d) )*((1-exp(-c*(x-t1)))/ (1-exp(-c*(t2-t1)))))^(1/d)" schnute

}
  \item{start}{
start values of the iteration process  
}

\item{mixed}{
FALSE/defalt for fixed model or TRUE for mixed model
}

\item{digits}{
number of digits in results (default=6)
}

\item{alpha}{
significant level of the confident intervals for parameters in the models
}

}

\value{
Returns coefficients of the models, t test for coefficients, R squared, adjusted R squared, AIC, BIC, and residuals of the model
}
\references{
KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.

TERRANCE J. QUINN II and RICHARD B. DERISO. Quantitative Fish Dynamics,
New York, Oxford, Oxford University Press, 1999. 
}
\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}




\seealso{
nls, nls2
}
\examples{
# weights of an Angus cow at ages from 8 to 108 months (Kaps and Lamberson, 2009)

weight=c(280,340,430,480,550,580,590,600,590,600)
age=c(8,12,24,36,48,60,72,84,96,108)

data1=data.frame(age, weight)

# linear
er1(data1, model=1)

# quadratic
er1(data1, model=2)

# linear plateau
er1(data1, model=3)

# quadratic plateau
er1(data1, model=4)

# two linear
er1(data1, model=5, start=c(250,6,2,50))

# exponential
er1(data1, model=6, start=c(250,0.05))

# logistic
er1(data1, model=7, start=c(600,4,0.05))

# van bertalanffy
er1(data1, model=8, start=c(600,2,0.05))

# brody
er1(data1, model=9, start=c(600,4,0.05))

# gompertz
er1(data1, model=10, start=c(600,4,0.05))

# richards
er1(data1, model=16, start=c(600,2,0.05,1.4))

# allometric
er1(data1, model=14)

# cubic
er1(data1, model=15)



# growth of Zagorje turkeys (Kaps and Lamberson, 2009)


weight=c(44,66,100,150,265,370,455,605,770)
age=c(1,7,14,21,28,35,42,49,56)

data2=data.frame(age,weight)

# two linear
er1(data2, model=5, start=c(25,6,10,20))

# gain weight measurements of turkey poults (Kaps and Lamberson, 2009)

methionine=c(80,85,90,95,100,105,110,115,120)
gain=c(102,115,125,133,140,141,142,140,142)

data3=data.frame(methionine, gain)

# linear
er1(data3, model=1)

# quadratic
er1(data3, model=2)

# linear plateau
er1(data3, model=3)

# quadratic plateau
er1(data3, model=4)

# lactation curve
 milk=c(25,24,26,28,30,31,27,26,25,24,23,24,22,21,22,
20,21,19,18,17,18,18,16,17,15,16,14)
 days=c(15,15,15,75,75,75,135,135,135,195,
195,195,255,255,255,315,315,315,375,375,375,435,435,435,495,495,495)
    
 data4=data.frame(days,milk)
	

er1(data4, model=11, start=c(16,0.25,0.004))

# ruminal degradation 
time=c(2,6,9,24,48,72,96)
deg=c(20,33,46,55,66,72,76)

data5=data.frame(time,deg)

er1(data5, model=12)

# logistic bi-compartmental (gas production)
time=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
gas=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    
data6=data.frame(time,gas)

er1(data6, model=13, start=c(19,4,0.025,0.004,5))

# Schnute model
#pacific halibut weight-age data of females (Terrance and Richard, 1999)
 
age=c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
19,20,21,22,23,24,28)
weight=c(1.7,2,3.9, 4.2,6.4,7.6,10.9,14.9,18.2,21.6,
25.4,28.8,30.9,	35.6,37.9,34.7,44.8,52.6,49.1,56.7,58.6,54.1)

halibut=data.frame(age,weight)


t1=min(halibut[,2])
t2=max(halibut[,2])

er1(halibut,model=17, start=c(a=t1,b=t2,c=0.15,d=-0.50))


}

