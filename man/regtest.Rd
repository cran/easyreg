\name{regtest}
\alias{regtest}

\title{
Test of models and parameters 
}
\description{
This function performs test of models and parameters
}
\usage{
regtest(data, model = 1, start = c(a = 1, b = 1, c = 1, d = 1,e = 1))
}

\arguments{
  \item{data}{
data is a data.frame
The first column contain explanatory variable, second column contain treatments and 
the third column contain the response variable
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
start values of iterations 
}  


}

\value{
Returns coefficients of the models, test for coefficients, AIC and BIC.
}

\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}



\seealso{
lm, ea1(easyanova package), pr2, regplot
}
\examples{
x=c(1,1,1,2,2,2,3,3,3,4,4,4)
y=c(5,5.3,6,8,8.9,12,14,18,25,25,29,32)
t=c("t1","t2","t3","t1","t2","t3","t1","t2","t3","t1","t2","t3")
data=data.frame(x,t,y)
# linear
regtest(data, model=1)
# quadratic
regtest(data, model=2)
# exponential
regtest(data, model=6)
# ... etc

}

