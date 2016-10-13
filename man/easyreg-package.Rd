\name{easyreg-package}
\alias{easyreg-package}
\alias{easyreg}
\docType{package}
\title{
Easy Regression
}
\description{
Performs analysis of regression in simple designs with quantitative treatments, including mixed models ans non linear models
}

\details{
\tabular{ll}{
Package: \tab easyreg\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-12\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Emmanuel Arnhold <emmanuelarnhold@yahoo.com.br>
}

\references{

KAPS, M. and LAMBERSON, W. R. Biostatistics for Animal Science: an introductory text. 2nd Edition. CABI Publishing, Wallingford, Oxfordshire, UK, 2009. 504p.

SAMPAIO, I. B. M. Estatistica aplicada a experimentacao animal. 3nd Edition. Belo Horizonte: Editora FEPMVZ, Fundacao de Ensino e Pesquisa em Medicina Veterinaria e Zootecnia, 2010. 264p.
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

# the growth of Zagorje turkeys (Kaps and Lamberson, 2009)

weight=c(44,66,100,150,265,370,455,605,770)
age=c(1,7,14,21,28,35,42,49,56)

data2=data.frame(age,weight)

# two linear
regplot(data2, model=5, start=c(25,6,10,20))

regplot(data2, model=5, start=c(25,6,10,20), digits=2)

# in other function
bl(data2)


}
