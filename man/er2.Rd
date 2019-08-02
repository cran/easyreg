\name{er2}
\alias{er2}

\title{
Analysis of polynomial regression
}
\description{
The function performs analysis of polynomial regression in simple designs with quantitative treatments. 
}
\usage{
er2(data, design = 1, list = FALSE, type = 2)
}

\arguments{
  \item{data}{
data is a data.frame

data frame with two columns, treatments and response (completely randomized design)

data frame with three columns, treatments, blocks and response (randomized block design)

data frame with four columns, treatments, rows, cols and response (latin square design)

data frame with five columns, treatments, square, rows, cols and response (several latin squares)
}
  \item{design}{

1 = completely randomized design

2 = randomized block design

3 = latin square design

4 = several latin squares

}
  \item{list}{
FALSE = a single response variable

TRUE = multivariable response
}
  \item{type}{
type is form of obtain sum of squares

1 = a sequential sum of squares 

2 = a partial sum of squares
}
}
\details{
The response  and the treatments must be numeric. Other variables can be numeric or factors.
}
\value{
Returns analysis of variance, models, t test for coefficients and R squared and adjusted R squared.
}
\references{
KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.

SAMPAIO, I. B. M. Estatistica aplicada a experimentacao animal. 3nd Edition. Belo Horizonte: Editora FEPMVZ, Fundacao de Ensino e Pesquisa em Medicina Veterinaria e Zootecnia, 2010. 264p.
}
\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}



\seealso{
lm, lme(package nlme), ea1(package easyanova), er1
}
\examples{
# analysis in completely randomized design
data(data1)
r1=er2(data1)
names(r1)
r1
r1[1]

# analysis in randomized block design
data(data2)
r2=er2(data2, design=2)
r2

# analysis in latin square design
data(data3)
r3=er2(data3, design=3)
r3

# analysis in several latin squares
data(data4)
r4=er2(data4, design=4)
r4

# data
treatments=rep(c(0.5,1,1.5,2,2.5,3), c(3,3,3,3,3,3))
r1=rnorm(18,60,3)
r2=r1*1:18
r3=r1*18:1
r4=r1*c(c(1:10),10,10,10,10,10,10,10,10)
data6=data.frame(treatments,r1,r2,r3, r4)

# use the argument list = TRUE
er2(data6, design=1, list=TRUE)

}

