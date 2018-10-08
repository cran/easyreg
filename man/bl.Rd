\name{bl}
\alias{bl}

\title{
Analysis of broken line regression
}
\description{
The function performs analysis of broken line regression 
}
\usage{
bl(data, model=1, alpha=0.05, xlab = "Explanatory Variable", ylab = "Response Variable", 
    position = 1, digits = 6, mean = TRUE, sd=FALSE, legend = TRUE, lty=2, 
col="dark blue", pch=20, xlim="default.x",ylim="default.y", ...)
}

\arguments{
  \item{data}{
data is a data.frame
The first column contain the treatments (explanatory variable) and the 
second column the response variable
}

 \item{model}{
model for analysis: 1=two linear; 2=linear plateau; 3=quadratic plateau; 4=model 1 with blocks random; 5=model 2 with effect of blocks; 6=model 3 with effect of blocks  
}

\item{alpha}{
significant level for cofidence intervals (parameters estimated)
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
number of digits (default=6)
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


#linear and quadratic plateau
x=c(0,1,2,3,4,5,6)
y=c(1,2,3,6.1,5.9,6,6.1)

data=data.frame(x,y)

bl(data,model=2, lty=1, col=1, digits=2, position=8)


# effect os blocks
x=c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)
y=c(4,12,9,20,16,25,21,31,28,42,33,46,33,46,34,44)
blocks=rep(c(1,2),8)

dat=data.frame(x,blocks,y)

bl(dat, 4)

bl(dat,5, sd=TRUE)

bl(dat,6, mean=FALSE)


}

