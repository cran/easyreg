\name{regplot}
\alias{regplot}

\title{
Plot data and equation
}
\description{
The function plot data and equation
}
\usage{
regplot(data, model=1, start=c(a=1,b=1,c=1,d=1,e=1), xlab="Explanatory Variable", 
ylab="Response Variable", position=1, digits=6, mean=TRUE, sd=FALSE, 
legend = TRUE, lty=2, col="dark blue", pch=20,  xlim="defalt.x",ylim="defalt.y",...)
}

\arguments{
  \item{data}{
data is a data.frame
The first column contain the treatments (explanatory variable) and the 
remaining column the response variable
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

}
  \item{start}{
start (iterations) values of model
}
  \item{xlab}{
names of variable x
}
  \item{ylab}{
names of variable y
}
  \item{position}{
position of equation in the graph

top=1

bottomright=2

bottom=3

bottomleft=4

left=5

topleft=6 (default)

topright=7

right=8

center=9
}

\item{digits}{
number of digits (defalt=6)
}

\item{mean}{
mean=TRUE (plot mean of data)
mean=FALSE (plot all data)
}

\item{sd}{
sd=FALSE (plot without standard deviation)
sd=TRUE (plot with standard deviation)
}

\item{legend}{
legend=TRUE (plot legend)
legend=FALSE (not plot legend)
}

\item{lty}{
line type
}

\item{col}{
line color
}

\item{pch}{
point type
}

\item{xlim}{
limits for x
}

\item{ylim}{
limits for y
}

\item{...}{
others graphical parameters (see par)
}

}


\references{
KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 
2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.
}
\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}



\seealso{
nls,er1,er2,bl
}
\examples{
# weights of Angus cow at ages from 8 to 108 months (Kaps and Lamberson, 2009)

weight=c(280,340,430,480,550,580,590,600,590,600)
age=c(8,12,24,36,48,60,72,84,96,108)

data1=data.frame(age, weight)

# linear
regplot(data1, model=1, digits=3, position=3, ylab="weight", xlab="age")

# quadratic
regplot(data1, model=2, digits=3, position=3, col=1, ylim=c(200,700))

# linear plateau
regplot(data1, model=3,ylab="weight", xlab="age", lty=5, col="dark green",
position=3, ylim=c(200,700), xlim=c(0,150), lwd=2)

# quadratic plateau
regplot(data1, model=4,ylab="weight", xlab="age")

# two linear
regplot(data1, model=5, start=c(250,6,2,50),digits=3, position=3 )

# exponential
regplot(data1, model=6, start=c(250,0.05))

# logistic
regplot(data1, model=7, start=c(600,4,0.05))

# van bertalanffy
regplot(data1, model=8, start=c(600,2,0.05))

# brody
regplot(data1, model=9, start=c(600,4,0.05))

# gompertz
regplot(data1, model=10, start=c(600,4,0.05))



# growth of Zagorje turkeys (Kaps and Lamberson, 2009)


weight=c(44,66,100,150,265,370,455,605,770)
age=c(1,7,14,21,28,35,42,49,56)

data2=data.frame(age,weight)

# two linear
regplot(data2, model=5, start=c(25,6,10,20))


# weight gain measurements of turkey poults (Kaps and Lamberson, 2009)

methionine=c(80,85,90,95,100,105,110,115,120)
gain=c(102,115,125,133,140,141,142,140,142)

data3=data.frame(methionine, gain)

# linear
regplot(data3, model=1)

# quadratic
regplot(data3, model=2)

# linear plateau
regplot(data3, model=3)

# quadratic plateau
regplot(data3, model=4)

# lactation curve
 milk=c(25,24,26,28,30,31,27,26,25,24,23,24,22,21,22,20,21,19,
18,17,18,18,16,17,15,16,14)
 days=c(15,15,15,75,75,75,135,135,135,195,195,195,255,255,255,315,
315,315,375,375,375,435,435,435,495,495,495)
    
data4=data.frame(days,milk)
	

regplot(data4, model=11, start=c(16,0.25,0.004))

# ruminal degradation 
time=c(2,6,9,24,48,72,96)
deg=c(20,33,46,55,66,72,76)

data5=data.frame(time,deg)

regplot(data5, model=12)

# logistic bi-compartmental (gas production)
time=c(0,12,24,36,48,60,72,84,96,108,120,144,168,192)
gas=c(0.002,3.8,8,14.5,16,16.5,17,17.4,17.9,18.1,18.8,19,19.2,19.3)
    
data6=data.frame(time,gas)

regplot(data6, model=13, start=c(19,4,0.025,0.004,5))

# multiple curves
time=c(0,12,24,48,64,72,96)
t1=c(36,48,59,72,85,86,87)
t2=c(14,25,36,49,59,65,72)
t3=c(55,78,86,87,86,87,88)

data=data.frame(time,t1,t2,t3)

regplot(data, model=12)
regplot(data, model=4)

# include standard deviation in graph
data(data1)

regplot(data1, sd=TRUE)


}

