\name{test_function}
\alias{test_function}
\title{Test functions for DetR}
\description{
     Functions to test the \code{cpp} codes in the package.
}
\usage{
  test_function()
}
\details{
	This is a series of \code{R} functions that, together, implement the \code{c++} codes used 
	in the package and which can be used to test those.
}
\examples{
n<-100
p<-5
set.seed(123) #for repoducibility.
Z<-matrix(rnorm(n*(p+1)),nc=p+1)
x<-Z[,1:p]
y<-Z[,p+1]
h<-ceiling((n+p)*0.6)
datao<-cbind(x,y)
test_R_0<-DetR:::test_fxOGK(x0=x,y0=y,cent_est='scaleTau2_test',scal_est='scaleTau2_test')
test_cpp<-DetR:::fxOGK(Data=datao,scale_est="scaleTau2",intercept=1,h=h,doCsteps=1)
####should be the same
sort(test_cpp$bestRaw)
sort(test_R_0$bestRaw)	
#############
test_R_1<-DetR:::test_Cstep(x=x,y=y,h=h,z0=test_R_0$bestRaw)
####should be the same
sort(test_R_1$bestCStep)
sort(test_cpp$bestCStep[1:h])
}
\author{
  Vakili Kaveh.
}
\keyword{multivariate}
\keyword{robust}
\keyword{deterministic}
