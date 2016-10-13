\name{bl}
\alias{bl}

\title{
Analysis of broken line regression
}
\description{
The function performs analysis of broken line regression 
}
\usage{
bl(data, xlab="Explanatory Variable", ylab="Response Variable", position=1,
 digits=6, mean=TRUE, legend=TRUE)
}

\arguments{
  \item{data}{
data is a data.frame
The first column contain the treatments (explanatory variable) and the 
second column the response variable
}

  \item{xlab}{
name of explanatory variable
}

  \item{ylab}{
name of response variable
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

\item{legend}{
legend=TRUE (plot legend)
legend=FALSE (not plot legend)
}

}

\value{
Returns coefficients of the models, t test for coefficients, knot (break point), R squared, adjusted R squared, AIC, BIC, residuals and shapiro-wilk test for residuals.
}

\references{
KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.
}


\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}



\seealso{
lm, ea1(easyanova package), er1
}
\examples{

# the growth of Zagorje turkeys (Kaps and Lamberson, 2009)

weight=c(44,66,100,150,265,370,455,605)
age=c(1,7,14,21,28,35,42,49)

data2=data.frame(age,weight)

# two linear
regplot(data2, model=5, start=c(25,6,10,20))

bl(data2, digits=2)



}

